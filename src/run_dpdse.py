# Here import the network details, create dpdse object instance and run dse
import cimpy
import numpy as np
import os
from pyvolt import network
from pyvolt import measurement
from pyvolt import nv_powerflow
from dpdse import DpDse, DpDse_Model, Line_Type
import config

xml_path = os.path.dirname(os.path.realpath(__file__)) + "/../data/"
xml_files = [os.path.join(xml_path, "10nodes_test_grid_v2.xml")]

# Read cim files and create new network.System object
res = cimpy.cim_import(xml_files, "cgmes_v2_4_15")
system = network.System()

# load cim data
base_apparent_power = 1  # MVA
Vbase = 12.66 # line-line voltage TODO: how to specify this uniformly? It is also used in dpdse class inside initialize_dse
system.load_cim_data(res['topology'], base_apparent_power)

# correct the system details, (if needed)
system.nodes[0].type = network.BusType.SLACK
system.nodes[0].voltage_pu = complex(1.0, 0)
system.nodes[0].voltage = complex(1.0*Vbase, 0)
system.nodes[4].type = network.BusType.PV
system.nodes[4].voltage_pu = complex(1.0, 0)
system.nodes[4].voltage = complex(1.0*Vbase, 0)
system.nodes[6].type = network.BusType.PV
system.nodes[6].voltage_pu = complex(1.0, 0)
system.nodes[6].voltage = complex(1.0*Vbase, 0)
system.nodes[8].type = network.BusType.PV
system.nodes[8].voltage_pu = complex(1.0, 0)
system.nodes[8].voltage = complex(1.0*Vbase, 0)

print('------------------Print elements UUIDs: --------------------')
for n in system.nodes:
    print(n.name, n.index, n.uuid)

for b in system.branches:
    print(b.start_node.index, b.end_node.index, b.uuid)

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
for node in results_pf.nodes:
    measurements_set.create_measurement(node.topology_node, measurement.ElemType.Node, measurement.MeasType.Vpmu_mag,
                                        np.absolute(node.voltage_pu), Pmu_mag_unc)
    measurements_set.create_measurement(node.topology_node, measurement.ElemType.Node, measurement.MeasType.Vpmu_phase,
                                        np.angle(node.voltage_pu), Pmu_phase_unc)
    measurements_set.create_measurement(node.topology_node, measurement.ElemType.Node, measurement.MeasType.V_mag,
                                        node.voltage_pu.real, V_unc)
    measurements_set.create_measurement(node.topology_node, measurement.ElemType.Node, measurement.MeasType.Sinj_real,
                                        node.power_pu.real, Sinj_unc)
    measurements_set.create_measurement(node.topology_node, measurement.ElemType.Node, measurement.MeasType.Sinj_imag,
                                        node.power_pu.imag, Sinj_unc)

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
                                        branch.power_pu.imag, S_unc)

measurements_set.meas_creation()


######################################## create an instance for dpdse and its config ####################################

dse_config = config.DSEConfig(config.DpDse_Mode.ONLINE, config.DpDse_Output.STORE, 10.00, config.DpDse_Smoother.No)
print(f"Mode: {dse_config.mode}, Output: {dse_config.output}, duration: {dse_config.duration}, smoother: {dse_config.smoother}")

gen_uuid = [gen_node.name for gen_node in system.get_ES_nodes()]
load_uuid = [load_node.name for load_node in system.get_EC_nodes()]

print("gen_nodes: ----: ", gen_uuid)

print("load_nodes: ----: ", load_uuid)

run_dpdse = DpDse(system, measurements_set, 1e-3, DpDse_Model.STANDARD, Line_Type.RL)
run_dpdse.initialize_dse()
run_dpdse.extract_u()


############################################## Run the DSE #############################################################
#TODO: create dse config object
if dse_config.mode == config.DpDse_Mode.ONLINE:
    # obtain measurements online and perform DSE for certain duration
    print("Under Construction!")
else:
    # obtain measurements from the file stored and perform DSE for certain duration
    print("Under Construction!")
