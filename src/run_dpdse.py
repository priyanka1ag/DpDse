# Here import the network details, create dpdse object instance and run dse
import cimpy
import numpy as np
import os
from pyvolt import network
from pyvolt import measurement
from pyvolt import nv_powerflow
from pyvolt import nv_state_estimator
from dpdse import DpDse, DpDse_Model, Line_Type
#from dpdse_ui import DpDse, DpDse_Model, Line_Type
import config



xml_path = os.path.dirname(os.path.realpath(__file__)) + "/../data/"
xml_files = [os.path.join(xml_path, "10nodes_test_grid.xml")]

# Read cim files and create new network.System object
res = cimpy.cim_import(xml_files, "cgmes_v2_4_15")
system = network.System()

# load cim data
base_apparent_power = 1  # MVA
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
c_node_index = [1, 2, 3, 4, 5, 6, 7, 8, 9] 
critical_nodes = [item for item in system.get_EC_nodes() if item.index in c_node_index]

# pass only required control inputs (gen voltage and critical load injection)
print("-----U Measurements---")
for node in results_pf.nodes:
    if node.topology_node.type == network.BusType.PV or node.topology_node.type == network.BusType.SLACK:
        print(f"Voltage: node uuid: {node.topology_node.uuid}, mag: , {np.absolute(node.voltage)}, ang: {np.angle(node.voltage)}")
        measurements_set.create_measurement(node.topology_node, measurement.ElemType.Node, measurement.MeasType.Vpmu_mag,
                                            np.absolute(node.voltage_pu), Pmu_mag_unc)
        measurements_set.create_measurement(node.topology_node, measurement.ElemType.Node, measurement.MeasType.Vpmu_phase,
                                            np.angle(node.voltage_pu), Pmu_phase_unc)
    
    elif node.topology_node.type == network.BusType.PQ and node.topology_node in critical_nodes:
        measurements_set.create_measurement(node.topology_node, measurement.ElemType.Node, measurement.MeasType.Ipmu_inj_mag,
                                            np.absolute(node.current_pu), Pmu_mag_unc)
        measurements_set.create_measurement(node.topology_node, measurement.ElemType.Node, measurement.MeasType.Ipmu_inj_phase,
                                            np.angle(node.current_pu), Pmu_phase_unc)
        print(f"Current: node uuid: {node.topology_node.uuid}, curr mag: {np.absolute(node.current)}, ang: {np.angle(node.current)}")

# following measurements for z
br_meas = [] # PyVolt branch object doesnt have index sadly!
load_vol_meas = [item for item in system.get_EC_nodes() if item.index in [4,6,8]] 
i = 0
print("-----Z Measurements---")
for br in results_pf.branches:
        if i in br_meas:
            print(f"Current: node uuid: {br.topology_branch.uuid}, curr mag: {np.absolute(br.current)}, ang: {np.angle(br.current)}, cmplx: {br.current}")
            measurements_set.create_measurement(br.topology_branch, measurement.ElemType.Branch, measurement.MeasType.Ipmu_mag ,
                                                np.absolute(br.current_pu), Pmu_mag_unc)
            measurements_set.create_measurement(br.topology_branch, measurement.ElemType.Branch, measurement.MeasType.Ipmu_phase,
                                                np.angle(br.current_pu), Pmu_phase_unc)
        i += 1

for node in results_pf.nodes:        
    if node.topology_node.type == network.BusType.PQ and node.topology_node in load_vol_meas:
        print(f"Voltage: node uuid: {node.topology_node.uuid}, mag: , {np.absolute(node.voltage)}, ang: {np.angle(node.voltage)}, cmplx: {node.voltage}")
        measurements_set.create_measurement(node.topology_node, measurement.ElemType.Node, measurement.MeasType.Vpmu_mag,
                                            np.absolute(node.voltage_pu), Pmu_mag_unc)
        measurements_set.create_measurement(node.topology_node, measurement.ElemType.Node, measurement.MeasType.Vpmu_phase,
                                            np.angle(node.voltage_pu), Pmu_phase_unc)
measurements_set.meas_creation()

######################################## create an instance for dpdse and its config ####################################

dse_config = config.DSEConfig(config.DpDse_Mode.OFFLINE, config.DpDse_Output.STORE, 10.00, config.DpDse_Smoother.No)
print(f"Mode: {dse_config.mode}, Output: {dse_config.output}, duration: {dse_config.duration}, smoother: {dse_config.smoother}")


#run_dpdse = DpDse(system, measurements_set, critical_nodes, 1e-3, DpDse_Model.STANDARD, Line_Type.RL)
run_dpdse = DpDse(system, measurements_set, 1e-3, DpDse_Model.STANDARD, Line_Type.RL)
run_dpdse.initialize_dse()
#run_dpdse.check_ss_consistency()



############################################## Run the DSE #############################################################
#TODO: create dse config object
if dse_config.mode == config.DpDse_Mode.ONLINE:
    # obtain measurements online and perform DSE for certain duration
    print("Under Construction!")
    run_dpdse.predict()
    #ONLINE: Check if new measurement received
    run_dpdse.correct()
    # stream results 
elif dse_config.mode == config.DpDse_Mode.OFFLINE:
    # obtain measurements from the file stored (time-series) and perform DSE for certain duration
    run_dpdse.predict()
    # some method to read and update next measurement line
    run_dpdse.correct()
    # store results as time-series
    print("Under Construction!")
