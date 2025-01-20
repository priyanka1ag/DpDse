# Here import the network details, create dpdse object instance and run dse
import cimpy
import numpy as np
import os
from pyvolt import network
from pyvolt import measurement
from pyvolt import nv_powerflow
from pyvolt import nv_state_estimator
from dpdse import DpDse, DpDse_Model, Line_Type
import config

xml_path = os.path.dirname(os.path.realpath(__file__)) + "/../data/"
xml_files = [os.path.join(xml_path, "seguro_net2.xml")]

# Read cim files and create new network.System object
res = cimpy.cim_import(xml_files, "cgmes_v2_4_15")
system = network.System()

# load cim data
base_apparent_power = 1  # MVA
Vbase = 10 # line-line voltage TODO: how to specify this uniformly? It is also used in dpdse class inside initialize_dse
system.load_cim_data(res['topology'], base_apparent_power)

# Check if voltage of Slack bus is 1+j0
for n in system.nodes:
    if n.type == network.BusType.SLACK:
        if n.voltage_pu == complex(0,0):
            n.voltage_pu = complex(1.0,0)
            n.voltage = complex(1.0*Vbase, 0)


print('------------------Print elements UUIDs: --------------------')
for n in system.nodes:
    print(f"Name: {n.name}, Index: {n.index}, UUID: {n.uuid}, Type: {n.type},  Voltage: {n.voltage_pu}")

for b in system.branches:
    print(f"Start node index: {b.start_node.index}, End node index: {b.end_node.index}, UUID: {b.uuid}, Resistance: {b.r}, Inductance: {b.x}, Capacitance: {b.c}")

# Execute power flow analysis
results_pf, num_iter = nv_powerflow.solve(system)

for node in results_pf.nodes:
    print(f"Bus {node.topology_node.index}, uuid: {node.topology_node.uuid}, voltage: {node.voltage}")


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



'''
for node in results_pf.nodes:
    if node.topology_node.type == network.BusType.PV or node.topology_node.type == network.BusType.SLACK:
        measurements_set.create_measurement(node.topology_node, measurement.ElemType.Node, measurement.MeasType.Vpmu_mag,
                                            np.absolute(node.voltage_pu), Pmu_mag_unc)
        measurements_set.create_measurement(node.topology_node, measurement.ElemType.Node, measurement.MeasType.Vpmu_phase,
                                            np.angle(node.voltage_pu), Pmu_phase_unc)
    elif node.topology_node.type == network.BusType.PQ:
        measurements_set.create_measurement(node.topology_node, measurement.ElemType.Node, measurement.MeasType.Ipmu_inj_mag,
                                            np.absolute(node.current_pu), Pmu_mag_unc)
        measurements_set.create_measurement(node.topology_node, measurement.ElemType.Node, measurement.MeasType.Ipmu_inj_phase,
                                            np.angle(node.current_pu), Pmu_phase_unc)

for branch in results_pf.branches:
    measurements_set.create_measurement(branch.topology_branch, measurement.ElemType.Branch, measurement.MeasType.Ipmu_mag,
                                        np.absolute(branch.current_pu), Pmu_mag_unc)
    measurements_set.create_measurement(branch.topology_branch, measurement.ElemType.Branch, measurement.MeasType.Ipmu_phase,
                                        np.angle(branch.current_pu), Pmu_phase_unc)
    measurements_set.create_measurement(branch.topology_branch, measurement.ElemType.Branch, measurement.MeasType.I_mag,
                                        branch.current_pu.real, I_unc)
    measurements_set.create_measurement(branch.topology_branch, measurement.ElemType.Branch, measurement.MeasType.S1_real,
                                        branch.power_pu.real, S_unc)
    measurements_set.create_measurement(branch.topology_branch, measurement.ElemType.Branch, measurement.MeasType.S1_imag,

'''
for node in results_pf.nodes:
    measurements_set.create_measurement(node.topology_node, measurement.ElemType.Node, measurement.MeasType.V_mag,
                                        np.absolute(node.voltage_pu), Pmu_mag_unc)
    
    #measurements_set.create_measurement(node.topology_node, measurement.ElemType.Node, measurement.MeasType.Vpmu_phase,
    #                                    np.angle(node.voltage_pu), Pmu_phase_unc)
    
    measurements_set.create_measurement(node.topology_node, measurement.ElemType.Node, measurement.MeasType.Sinj_real,
                                        np.real(node.power_pu), Pmu_mag_unc)
    
    measurements_set.create_measurement(node.topology_node, measurement.ElemType.Node, measurement.MeasType.Sinj_imag,
                                        np.imag(node.power_pu), Pmu_mag_unc)

measurements_set.meas_creation()


# Perform state estimation
state_estimation_results = nv_state_estimator.DsseCall(system, measurements_set)

# Print node voltages
print("state_estimation_results.voltages: ")
for node in state_estimation_results.nodes:
    print('{}={}'.format(node.topology_node.uuid, node.voltage))




######################################## create an instance for dpdse and its config ####################################

dse_config = config.DSEConfig(config.DpDse_Mode.ONLINE, config.DpDse_Output.STORE, 10.00, config.DpDse_Smoother.No)
print(f"Mode: {dse_config.mode}, Output: {dse_config.output}, duration: {dse_config.duration}, smoother: {dse_config.smoother}")

gen_uuid = [gen_node.name for gen_node in system.get_ES_nodes()]
load_uuid = [load_node.name for load_node in system.get_EC_nodes()]

run_dpdse = DpDse(system, measurements_set, 1e-3, DpDse_Model.STANDARD, Line_Type.RL)
run_dpdse.initialize_dse()
run_dpdse.check_ss_consistency()


############################################## Run the DSE #############################################################
#TODO: create dse config object
if dse_config.mode == config.DpDse_Mode.ONLINE:
    # obtain measurements online and perform DSE for certain duration
    print("Under Construction!")
else:
    # obtain measurements from the file stored and perform DSE for certain duration
    print("Under Construction!")


