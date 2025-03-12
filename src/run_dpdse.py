# Here import the network details, create dpdse object instance and run dse
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
import mqtt_fetch 
import csv

from seguro.common.store import Client, Event
from seguro.common import store, job


xml_path = os.path.dirname(os.path.realpath(__file__)) + "/../data/"
xml_files = [os.path.join(xml_path, "10nodes_test_grid.xml")]

# Read cim files and create new network.System object
res = cimpy.cim_import(xml_files, "cgmes_v2_4_15")
system = network.System()

# load cim data
base_apparent_power = 100  # MVA
Vbase = 12.66  # line-line voltage TODO: how to specify this uniformly? It is also used in dpdse class inside initialize_dse
system.load_cim_data(res['topology'], base_apparent_power)

# correct the system details, (if needed)
system.nodes[0].type = network.BusType.SLACK
system.nodes[0].voltage_pu = complex(0, -0.9996) 
system.nodes[0].voltage = complex(0, -0.9996)*Vbase 

print('------------------Print elements UUIDs: --------------------')
for n in system.nodes:
    print(n.name, n.index, n.uuid, n.type)

for b in system.branches:
    print(b.start_node.index, b.end_node.index, b.uuid, b.base_current, b.r, b.x/(2*np.pi*50))

# Execute power flow analysis
results_pf, num_iter = nv_powerflow.solve(system)


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
c_node_index = [1, 2, 3, 4, 5, 6, 7, 8, 9] 
critical_nodes = [item for item in system.get_EC_nodes() if item.index in c_node_index]
pq_nodes = [1, 2, 3, 4, 5, 6, 7, 8, 9]
curr_nodes = []
# pass only required control inputs (gen voltage and critical load injection)
print("-----U Measurements---")
for node in results_pf.nodes:
    if node.topology_node.type == network.BusType.PV or node.topology_node.type == network.BusType.SLACK:
        print(f"Voltage: node uuid: {node.topology_node.uuid}, name: {node.topology_node.name}, index: {node.topology_node.index}, mag: , {np.absolute(node.voltage)}, ang: {np.angle(node.voltage)}, cmplx: {node.voltage}")
        measurements_set.create_measurement(node.topology_node, measurement.ElemType.Node, measurement.MeasType.Vpmu_mag,
                                            np.absolute(node.voltage_pu), Pmu_mag_unc)
        measurements_set.create_measurement(node.topology_node, measurement.ElemType.Node, measurement.MeasType.Vpmu_phase,
                                            np.angle(node.voltage_pu), Pmu_phase_unc)
    
    elif node.topology_node.type == network.BusType.PQ and node.topology_node in critical_nodes and node.topology_node.index in curr_nodes:
        measurements_set.create_measurement(node.topology_node, measurement.ElemType.Node, measurement.MeasType.Ipmu_inj_mag,
                                            np.absolute(node.current_pu), Pmu_mag_unc)
        measurements_set.create_measurement(node.topology_node, measurement.ElemType.Node, measurement.MeasType.Ipmu_inj_phase,
                                            np.angle(node.current_pu), Pmu_phase_unc)
        print(f"Current: node uuid: {node.topology_node.uuid}, curr mag: {np.absolute(node.current)}, ang: {np.angle(node.current)}, cmplx: {node.current}")
    elif node.topology_node.type == network.BusType.PQ and node.topology_node in critical_nodes and node.topology_node.index in pq_nodes:
        measurements_set.create_measurement(node.topology_node, measurement.ElemType.Node, measurement.MeasType.Sinj_real,
                                            np.real(node.power_pu), Pmu_mag_unc)
        measurements_set.create_measurement(node.topology_node, measurement.ElemType.Node, measurement.MeasType.Sinj_imag,
                                            np.imag(node.power_pu), Pmu_phase_unc)
        print(f"Power: node uuid: {node.topology_node.uuid}, P: {np.real(node.power)}, Q: {np.imag(node.power)}, i_inj_cmplx: {(node.current)}, vl_cmplx: {node.voltage}")



# following measurements for z
br_meas_pmu = [0, 4, 6] # PyVolt branch object doesnt have index sadly! 0, 4, 6
br_meas_scada = [2, 8, 5]
vol_meas = [1, 8, 5, 6, 7] # 4, 6, 8
vol_mag_scada = []
load_vol_meas = [item for item in system.get_EC_nodes() if item.index in vol_meas] 
i = 0
print("-----Z Measurements---")
for br in results_pf.branches:
        if i in br_meas_pmu:
            print(f"Current: node uuid: {br.topology_branch.uuid}, curr_pu_mag: {np.absolute(br.current_pu)} curr mag: {np.absolute(br.current)}, ang: {np.angle(br.current)}, cmplx: {br.current}")
            measurements_set.create_measurement(br.topology_branch, measurement.ElemType.Branch, measurement.MeasType.Ipmu_mag ,
                                                np.absolute(br.current_pu), Pmu_mag_unc)
            measurements_set.create_measurement(br.topology_branch, measurement.ElemType.Branch, measurement.MeasType.Ipmu_phase,
                                                np.angle(br.current_pu), Pmu_phase_unc)
        if i in br_meas_scada:
            print(f"Current: node uuid: {br.topology_branch.uuid}, curr_pu_mag: {np.absolute(br.current_pu)} curr mag: {np.absolute(br.current)}, ang: {np.angle(br.current)}, cmplx: {br.current}")
            measurements_set.create_measurement(br.topology_branch, measurement.ElemType.Branch, measurement.MeasType.I_mag ,
                                                np.absolute(br.current_pu), I_unc)
        i += 1

for node in results_pf.nodes:        
    if node.topology_node.type == network.BusType.PQ and node.topology_node in load_vol_meas:
        print(f"Voltage: node uuid: {node.topology_node.uuid}, name: {node.topology_node.name}, index: {node.topology_node.index}, mag: , {np.absolute(node.voltage)}, ang: {np.angle(node.voltage)}, cmplx: {node.voltage}")
        if node.topology_node.index in vol_mag_scada:
            measurements_set.create_measurement(node.topology_node, measurement.ElemType.Node, measurement.MeasType.V_mag,
                                                np.absolute(node.voltage_pu), V_unc)
        else:
            measurements_set.create_measurement(node.topology_node, measurement.ElemType.Node, measurement.MeasType.Vpmu_mag,
                                                np.absolute(node.voltage_pu), Pmu_mag_unc)
            measurements_set.create_measurement(node.topology_node, measurement.ElemType.Node, measurement.MeasType.Vpmu_phase,
                                                np.angle(node.voltage_pu), Pmu_phase_unc)
measurements_set.meas_creation()

######################################## create an instance for dpdse and its config ####################################

dse_config = config.DSEConfig(config.DpDse_Mode.OFFLINE, config.DpDse_Output.STORE, 10.00, config.DpDse_Smoother.No)
print(f"Mode: {dse_config.mode}, Output: {dse_config.output}, duration: {dse_config.duration}, smoother: {dse_config.smoother}")

run_dpdse = DpDse(system, measurements_set, 0.02, Line_Type.RL)
run_dpdse.initialize_dse()
map_u = {(m.element.uuid, m.meas_type): m for m in run_dpdse.get_meas_u().measurements}
map_z = {(m.element.uuid, m.meas_type): m for m in run_dpdse.get_meas_z().measurements}
#print("map_u keys: ", map_u.keys())
#print("map_z keys: ", map_z.keys())
#run_dpdse.check_ss_consistency()



############################################## Run the DSE #############################################################
#TODO: create dse config object
if dse_config.mode == config.DpDse_Mode.ONLINE:
    # obtain measurements online and perform DSE for certain duration
    print("Under Construction!")
    #ONLINE: Check if new measurement received
    topic = "network_10_nodes"
    mqtt_fetch.rerun_every_minute(topic, run_dpdse)
    # stream results 

elif dse_config.mode == config.DpDse_Mode.OFFLINE:
    # obtain measurements from the file stored (time-series) and perform DSE for certain duration
    # some method to read and update next measurement line
    meas_files = os.path.join(xml_path, "updated_data_Bus9_Fault.csv")
    df = pd.read_csv(meas_files)
    value_columns = df.columns[6:] # only time-series measurement values

    print("len of time-series: ", np.shape(value_columns))
    columns_names = []
    for key in run_dpdse.states_output.keys():
        columns_names.extend([f"{key}_Real", f"{key}_Imag", f"{key}_Mag", f"{key}_Phase"])

    est_values = []

    for t in value_columns:
        each_row = []
        print(f"--------------Updating matters to {t}-----------------")
        # Iterate over all columns except 'uid' and 'type'
        for _, row in df.iterrows():
            key = (str(row['UUID']), measurement.MeasType(int(row['Type'])))
            
            new_meas_val = row[t]
            if row['Unit'] == 'V' or row['Unit'] == 'A':
                new_meas_val = row[t]/1000
            if measurement.MeasType(int(row['Type'])) in [2, 3, 4, 5]: # converting power measurements to per phase values
                new_meas_val = row[t]/3
            if key in map_u:
                mea = map_u[key]
                run_dpdse.update_measurement(str(row['UUID']), measurement.MeasType(int(row['Type'])), new_meas_val, map_u, value_in_pu=row['Pu']) 
            elif key in map_z:
                mea = map_z[key]
                run_dpdse.update_measurement(str(row['UUID']), measurement.MeasType(int(row['Type'])), new_meas_val, map_z, value_in_pu=row['Pu']) 
            else:
                pass
                #print(f"Matter with uid={row['UUID']} and type={row['Type']} not found.")
        

        run_dpdse.predict()
        run_dpdse.correct()

        # Extract values from the dictionary and append to the row
        for key, val in run_dpdse.states_output.items():
            for array in val:
                each_row.append(array.item()) # Get values (real, imag, mag, phase)
        est_values.append(each_row)


    output_df = pd.DataFrame(est_values, columns=columns_names)
    output_df.to_csv('myDSEoutput.csv', index=False)

    # create seguro store client
    store = Client()
    
    # Put file into storage
    store.put_file("myDSEoutput.csv", "myDSEoutput.csv")
 