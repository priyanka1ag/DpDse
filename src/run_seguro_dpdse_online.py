# SPDX-FileCopyrightText: 2023-2024 Steffen Vogel, OPAL-RT Germany GmbH
# SPDX-License-Identifier: Apache-2.0

import sys
import time
import argparse
import environ
import logging
from functools import partial

import cimpy
import numpy as np
import os
from pyvolt import network
from pyvolt import measurement
from pyvolt import nv_powerflow
from pyvolt import nv_state_estimator
from dpdse import DpDse, Line_Type
import config
import pandas as pd
import csv
import re

from seguro.common import broker, config
from villas.node.sample import Sample

env = environ.Env()

# Authentication
TOPIC = env.str("TOPIC", "data/measurements/mp1")
TOPIC_PROCESSED = env.str(
    "TOPIC", "data/measurements/processed_by_streaming_worker/mp1"
)
RATE = env.float("RATE", 10.0)
VALUES = env.int("VALUES", 6)
BLOCK_INTERVAL = env.str("BLOCK_INTERVAL", "1m")


def new_samples(args, b: broker.Client, topic: str, samples: list[Sample]):
    for sample in samples:
        print(topic, sample)
        
        # Here update the measurements and run correct step.
        run_dpdse = args[0]
        map_u = args[1]
        map_z = args[2]
        
        for key_u, m_u in map_u.items():
            uuid = key_u[0]
            meas_type = key_u[1]
            value = m_u.meas_value_ideal # replacing with PF original value
            #print("key u: ", key_u, value)
            run_dpdse.update_measurement(str(uuid), meas_type, value, map_u, value_in_pu=True) 
        for key_z, m_z in map_z.items():
            uuid = key_z[0]
            meas_type = key_z[1]
            value = m_z.meas_value_ideal # replacing with PF original value
            #print("key z: ", key_z, value)
            run_dpdse.update_measurement(str(uuid), meas_type, value, map_z, value_in_pu=True) 
        
        
        run_dpdse.correct()
        print("CORRECT")

        for i, _ in enumerate(sample.data):
            sample.data[i] *= i



def main() -> int:
    xml_path = os.path.dirname(os.path.realpath(__file__)) + "/../data/"
    xml_files = [os.path.join(xml_path, "seguro_split_net2_DPSimWorking.xml")]

    # Read cim files and create new network.System object
    res = cimpy.cim_import(xml_files, "cgmes_v2_4_15")
    system = network.System()

    # load cim data
    base_apparent_power = 100  # MVA
    Vbase = 10 # line-line voltage TODO: how to specify this uniformly? It is also used in dpdse class inside initialize_dse
    system.load_cim_data(res['topology'], base_apparent_power)

    # Check if voltage of Slack bus is 1+j0
    for n in system.nodes:
        #    if n.voltage_pu == complex(0,0):
                n.voltage_pu = complex(1.0,0)
                n.voltage = complex(1.0*Vbase, 0)

    # Execute power flow analysis
    results_pf, num_iter = nv_powerflow.solve(system)


    print('------------------Print Power Flow : --------------------')
    print(f'------------------ Power Flow solved in : {num_iter} iterations --------------------')
    print(f'------------------ Power Flow NODES --------------------')
    for node in results_pf.nodes:
        print(f"Uuid: {node.topology_node.uuid}, name: {node.topology_node.name}, V: {node.voltage}, V_mag: {np.absolute(node.voltage)}, V_ang: {np.angle(node.voltage, deg = True)}, S: {node.power}")

    print(f'------------------ Power Flow BRANCHES --------------------')
    for br in results_pf.branches:
        print(f"uuid: {br.topology_branch.uuid}, curr: {br.current}, curr_pu: {br.current_pu}, curr_mag: {np.absolute(br.current)}, curr_ang: {np.angle(br.current, deg = True)} ")


    #################################### Step-2 Declaring information about measurement devices ####################################
    """ Write here the percent uncertainties of the measurements"""
    V_unc = 0
    I_unc = 0
    Sinj_unc = 0
    S_unc = 0
    Pmu_mag_unc = 0
    Pmu_phase_unc = 0

    # Create measurements data structures
    """first create measurement object for required measurements + control inputs"""

    measurements_set = measurement.MeasurementSet()

    # Create measurements data structures
    """first create measurement object for required measurements + control inputs"""

    measurements_set = measurement.MeasurementSet()
    c_node_index = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 20, 21] 
    critical_nodes = [item for item in system.get_EC_nodes() if item.index in c_node_index]
    pq_nodes = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 20, 21]
    curr_nodes = []
    # pass only required control inputs (gen voltage and critical load injection)
    #print("-----U Measurements---")

    for node in results_pf.nodes:
        if node.topology_node.type == network.BusType.PV or node.topology_node.type == network.BusType.SLACK:
            #print(f"Voltage: node uuid: {node.topology_node.uuid}, name: {node.topology_node.name}, index: {node.topology_node.index}, mag: , {np.absolute(node.voltage)}, ang: {np.angle(node.voltage)}, cmplx: {node.voltage}")
            measurements_set.create_measurement(node.topology_node, measurement.ElemType.Node, measurement.MeasType.Vpmu_mag,
                                                np.absolute(node.voltage_pu), Pmu_mag_unc)
            measurements_set.create_measurement(node.topology_node, measurement.ElemType.Node, measurement.MeasType.Vpmu_phase,
                                                np.angle(node.voltage_pu), Pmu_phase_unc)
        
        elif node.topology_node.type == network.BusType.PQ and node.topology_node in critical_nodes and node.topology_node.index in curr_nodes:
            measurements_set.create_measurement(node.topology_node, measurement.ElemType.Node, measurement.MeasType.Ipmu_inj_mag,
                                                np.absolute(node.current_pu), Pmu_mag_unc)
            measurements_set.create_measurement(node.topology_node, measurement.ElemType.Node, measurement.MeasType.Ipmu_inj_phase,
                                                np.angle(node.current_pu), Pmu_phase_unc)
            #print(f"Current: node uuid: {node.topology_node.uuid}, curr mag: {np.absolute(node.current)}, ang: {np.angle(node.current)}, cmplx: {node.current}")
        elif node.topology_node.type == network.BusType.PQ and node.topology_node in critical_nodes and node.topology_node.index in pq_nodes:
            measurements_set.create_measurement(node.topology_node, measurement.ElemType.Node, measurement.MeasType.Sinj_real,
                                                np.real(node.power_pu), Pmu_mag_unc)
            measurements_set.create_measurement(node.topology_node, measurement.ElemType.Node, measurement.MeasType.Sinj_imag,
                                                np.imag(node.power_pu), Pmu_phase_unc)
            #print(f"Power: node uuid: {node.topology_node.uuid}, P: {np.real(node.power)}, Q: {np.imag(node.power)}, i_inj_cmplx: {(node.current)}, vl_cmplx: {node.voltage}")



    # following measurements for z
    br_meas_pmu = [1, 3, 5, 8, 14, 16, 18, 19, 20] # PyVolt branch object doesnt have index sadly! 0, 4, 6
    br_meas_scada = []
    vol_meas = [1, 2, 8, 12, 14, 21] # 4, 6, 8
    vol_mag_scada = []
    load_vol_meas = [item for item in system.get_EC_nodes() if item.index in vol_meas] 
    i = 0
    #print("-----Z Measurements---")
    for br in results_pf.branches:
            if i in br_meas_pmu:
                #print(f"Current: node uuid: {br.topology_branch.uuid}, curr_pu_mag: {np.absolute(br.current_pu)} curr mag: {np.absolute(br.current)}, ang: {np.angle(br.current)}, cmplx: {br.current}")
                measurements_set.create_measurement(br.topology_branch, measurement.ElemType.Branch, measurement.MeasType.Ipmu_mag ,
                                                    np.absolute(br.current_pu), Pmu_mag_unc)
                measurements_set.create_measurement(br.topology_branch, measurement.ElemType.Branch, measurement.MeasType.Ipmu_phase,
                                                    np.angle(br.current_pu), Pmu_phase_unc)
            if i in br_meas_scada:
                #print(f"Current: node uuid: {br.topology_branch.uuid}, curr_pu_mag: {np.absolute(br.current_pu)} curr mag: {np.absolute(br.current)}, ang: {np.angle(br.current)}, cmplx: {br.current}")
                measurements_set.create_measurement(br.topology_branch, measurement.ElemType.Branch, measurement.MeasType.I_mag ,
                                                    np.absolute(br.current_pu), I_unc)
            i += 1

    for node in results_pf.nodes:        
        if node.topology_node.type == network.BusType.PQ and node.topology_node in load_vol_meas:
            #print(f"Voltage: node uuid: {node.topology_node.uuid}, name: {node.topology_node.name}, index: {node.topology_node.index}, mag: , {np.absolute(node.voltage)}, ang: {np.angle(node.voltage)}, cmplx: {node.voltage}")
            if node.topology_node.index in vol_mag_scada:
                measurements_set.create_measurement(node.topology_node, measurement.ElemType.Node, measurement.MeasType.V_mag,
                                                    np.absolute(node.voltage_pu), V_unc)
            else:
                measurements_set.create_measurement(node.topology_node, measurement.ElemType.Node, measurement.MeasType.Vpmu_mag,
                                                    np.absolute(node.voltage_pu), Pmu_mag_unc)
                measurements_set.create_measurement(node.topology_node, measurement.ElemType.Node, measurement.MeasType.Vpmu_phase,
                                                    np.angle(node.voltage_pu), Pmu_phase_unc)

    measurements_set.meas_creation()



    # Perform state estimation
    state_estimation_results = nv_state_estimator.DsseCall(system, measurements_set)

    # Print node voltages
    print("state_estimation_results.voltages: ")

    ######################################## create an instance for dpdse and its config ####################################

    gen_uuid = [gen_node.name for gen_node in system.get_ES_nodes()]
    load_uuid = [load_node.name for load_node in system.get_EC_nodes()]

    run_dpdse = DpDse(system, measurements_set, 0.0001, Line_Type.RL)
    run_dpdse.initialize_dse()
    map_u = {(m.element.uuid, m.meas_type): m for m in run_dpdse.get_meas_u().measurements}
    map_z = {(m.element.uuid, m.meas_type): m for m in run_dpdse.get_meas_z().measurements}
    
    args = [run_dpdse, map_u, map_z]

    b = broker.Client("example-streaming-worker")

    b.subscribe_samples(TOPIC, partial(new_samples, args))


    # run first predict (as otherwise, the program starts with first correct, and then it leads to error)
    run_dpdse.predict()

    while True:
        try:
            time.sleep(1)
            run_dpdse.predict()
            print("PREDICTED")
        except KeyboardInterrupt:
            break

    return 0


if __name__ == "__main__":
    sys.exit(main())