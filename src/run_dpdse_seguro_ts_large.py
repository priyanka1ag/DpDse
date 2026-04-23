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
import csv
import re
import time


from seguro.common.store import Client, Event
from seguro.common import store, job


#################################### Step-1 Read CIM file and run power flow ####################################
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


print(f'------------------ Power Flow solved in : {num_iter} iterations --------------------')

#################################### Step-2 Declare information about measurement devices ####################################
""" Write here the percent uncertainties of the measurements"""
V_unc = 0.2
I_unc = 0.2
Sinj_unc = 0.2
S_unc = 0.2
Pmu_mag_unc = 0.2
Pmu_phase_unc = 0.0001

# Create measurements data structures
"""first create measurement object for required measurements + control inputs"""

measurements_set = measurement.MeasurementSet()
vl_nodes = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 20, 21]
br_meas_pmu = [1] # PyVolt branch object doesnt have index sadly! 0, 4, 6
pq_nodes = []
# pass only required control inputs (gen voltage and load node voltages or power injections)

for node in results_pf.nodes:
    if node.topology_node.type == network.BusType.PV or node.topology_node.type == network.BusType.SLACK:
        #print(f"Voltage: node uuid: {node.topology_node.uuid}, name: {node.topology_node.name}, index: {node.topology_node.index}, mag: , {np.absolute(node.voltage)}, ang: {np.angle(node.voltage)}, cmplx: {node.voltage}")
        measurements_set.create_measurement(node.topology_node, measurement.ElemType.Node, measurement.MeasType.Vpmu_mag,
                                            np.absolute(node.voltage_pu), Pmu_mag_unc)
        measurements_set.create_measurement(node.topology_node, measurement.ElemType.Node, measurement.MeasType.Vpmu_phase,
                                            np.angle(node.voltage_pu), Pmu_phase_unc)
        #print(f"Current: node uuid: {node.topology_node.uuid}, curr mag: {np.absolute(node.current)}, ang: {np.angle(node.current)}, cmplx: {node.current}")
    elif node.topology_node.type == network.BusType.PQ and node.topology_node.index in pq_nodes:
        measurements_set.create_measurement(node.topology_node, measurement.ElemType.Node, measurement.MeasType.Sinj_real,
                                            np.real(node.power_pu), Pmu_mag_unc)
        measurements_set.create_measurement(node.topology_node, measurement.ElemType.Node, measurement.MeasType.Sinj_imag,
                                            np.imag(node.power_pu), Pmu_phase_unc)
        #print(f"Power: node uuid: {node.topology_node.uuid}, P: {np.real(node.power)}, Q: {np.imag(node.power)}, i_inj_cmplx: {(node.current)}, vl_cmplx: {node.voltage}")
    elif node.topology_node.type == network.BusType.PQ  and node.topology_node.index in vl_nodes:
        #print(f"Voltage: node uuid: {node.topology_node.uuid}, name: {node.topology_node.name}, index: {node.topology_node.index}, mag: , {np.absolute(node.voltage)}, ang: {np.angle(node.voltage)}, cmplx: {node.voltage}")
        measurements_set.create_measurement(node.topology_node, measurement.ElemType.Node, measurement.MeasType.Vpmu_mag,
                                            np.absolute(node.voltage_pu), Pmu_mag_unc)
        measurements_set.create_measurement(node.topology_node, measurement.ElemType.Node, measurement.MeasType.Vpmu_phase,
                                           np.angle(node.voltage_pu), Pmu_phase_unc)

i = 0
for br in results_pf.branches:
        if i in br_meas_pmu:
            #print(f"Current: node uuid: {br.topology_branch.uuid}, curr_pu_mag: {np.absolute(br.current_pu)} curr mag: {np.absolute(br.current)}, ang: {np.angle(br.current)}, cmplx: {br.current}")
            measurements_set.create_measurement(br.topology_branch, measurement.ElemType.Branch, measurement.MeasType.Ipmu_mag ,
                                                np.absolute(br.current_pu), Pmu_mag_unc)
            measurements_set.create_measurement(br.topology_branch, measurement.ElemType.Branch, measurement.MeasType.Ipmu_phase,
                                                np.angle(br.current_pu), Pmu_phase_unc)
        i += 1

measurements_set.meas_creation()

######################################## Step-3 create an instance for dpdse ####################################

gen_uuid = [gen_node.name for gen_node in system.get_ES_nodes()]
load_uuid = [load_node.name for load_node in system.get_EC_nodes()]

run_dpdse = DpDse(system, measurements_set, 0.01, Line_Type.RL)
run_dpdse.initialize_dse()
map_u = {(m.element.uuid, m.meas_type): m for m in run_dpdse.get_meas_u().measurements}
map_z = {(m.element.uuid, m.meas_type): m for m in run_dpdse.get_meas_z().measurements}

nodes_uid = [n.uuid for n in system.nodes]
branches_uid = [b.uuid for b in system.branches]

run_dpdse.check_ss_consistency()
Act = run_dpdse.get_Act()
eigenvalues = np.linalg.eigvals(Act)
frequencies = np.abs(np.real(eigenvalues))
highest_frequency = np.max(frequencies)
#print("highest frequency: ", highest_frequency)
print("time-step suggested: ", 1/(highest_frequency))




############################################## Run the DSE #############################################################

# obtain measurements from the file stored and perform DSE for certain duration
# obtain measurements from the file stored (time-series) and perform DSE for certain duration
# some method to read and update next measurement line
#meas_files = os.path.join(xml_path, "dynamic_SEGUROGrid_Fault_usable.csv")

# create seguro store client
#store = Client()

meas_file_name = "seguro_split_net2.csv"
#store.get_file(meas_file_name, "data/results/state-estimation/dynamic/seguro_split_net2") # TODO: path where the DPSim results are stored
meas_file_path = os.path.dirname(os.path.realpath(__file__)) + "/" + meas_file_name 
df_rect = pd.read_csv(meas_file_path)

# To change measurement from rect to polar
pattern = re.compile(r'(?P<element>.+?)\.(?P<type>[IV])\.(?P<part>re|im)')

# Step 3: Group real and imaginary columns
components = {}
for col in df_rect.columns:
    match = pattern.fullmatch(col)
    if match:
        elem = match.group('element')
        signal_type = match.group('type')
        part = match.group('part')
        key = f"{elem}.{signal_type}"
        if key not in components:
            components[key] = {}
        components[key][part] = col

# Step 4: Create list to collect new columns
new_columns = {}

# Step 5: Process each pair to compute magnitude and angle
for key, parts in components.items():
    if 're' in parts and 'im' in parts:
        re_col = parts['re']
        im_col = parts['im']

        complex_vals = df_rect[re_col].values + 1j * df_rect[im_col].values

        # New column names
        mag_col = re_col.replace('.re', '.mag')
        ang_col = im_col.replace('.im', '.ang')

        new_columns[mag_col] = np.abs(complex_vals)
        new_columns[ang_col] = np.angle(complex_vals)

# Step 6: Drop only used real/imag columns
re_im_cols_to_drop = [parts['re'] for parts in components.values() if 're' in parts and 'im' in parts] + \
                     [parts['im'] for parts in components.values() if 're' in parts and 'im' in parts]

df_cleaned = df_rect.drop(columns=re_im_cols_to_drop)

# Step 7: Combine everything with concat (avoids fragmentation warning)
new_df = pd.concat([df_cleaned, pd.DataFrame(new_columns)], axis=1)

# Step 8: Save to CSV - polar measurements
new_df.to_csv('dynamic_SEGUROGrid_usable.csv', index=False)



# Read measurements from polar file
value_columns = new_df.columns # only time-series measurement values

new_columns1 = {}

# suffix map to change the column name from the meas file received from DP-Sim
# Note that below, it does not cover all measurement types - ex Vmag, Imag, or S1_real, S1_imag, S2_real, S2_imag etc
suffix_map = {
('node', 'S.re'): 2,
('node', 'S.im'): 3,
('branch', 'S.re'): 4,
('branch', 'S.im'): 5,
('node', 'V.mag'): 7, 
('node', 'V.ang'): 8,
('node', 'I.mag'): 13,
('node', 'I.ang'): 14,
('branch', 'I.mag'): 9,
('branch', 'I.ang'): 10
# Add more if needed
}

for col in new_df.columns:
    parts = col.split('.')
    if len(parts) >= 3:
        uuid = parts[0].strip()
        suffix = '.'.join(parts[1:])  # e.g., 'S.re'
        if uuid in nodes_uid:
            key = ('node', suffix)
        elif uuid in branches_uid:
            key = ('branch', suffix)
        else:
            continue  # or log a warning
        mapped_number = suffix_map.get(key)
        if mapped_number is not None:
            new_columns1[col] = f"{uuid}.{mapped_number}"

# Rename columns
new_df.rename(columns=new_columns1, inplace=True)

columns_names = []
for key in run_dpdse.states_output.keys():
    columns_names.extend([f"{key}_re", f"{key}_im", f"{key}_Mag", f"{key}_Phase", f"{key}_re_var", f"{key}_im_var"])

n_rows = new_df.shape[0]
est_values = []

# record noisy measurements of voltage at tillmannshof
v_tillmannshof_mag = np.empty(n_rows) 
v_tillmannshof_ang = np.empty(n_rows)

# record WLS estimator output
v_tillmannshof_re_wls = np.empty(n_rows)
v_tillmannshof_im_wls = np.empty(n_rows)
i_ltg_l07H_t02h_re_wls = np.empty(n_rows)
i_ltg_l07H_t02h_im_wls = np.empty(n_rows)

i = 0

for _, row in new_df.iterrows():
#while i < 4:
    each_row = []   
    
    for col in new_df.columns:
        value = row[col]
        if pd.notnull(value):
            try:
                uuid, num = col.split('.')
                # convert the real/reactive powers to MVar/MW and voltages/currents to kV/kA and degrees to rads
                if int(num) in [2, 3, 4, 5, 11, 12]:
                    value = -value /1000000
                elif int(num) in [1, 6, 7, 9, 13]:
                    value = value /(1000)
                elif int(num) in [10]: # branch current directions are opposite - magnitude is same but phase is 180 shifted
                     value = value + np.pi

                key = (str(uuid), measurement.MeasType(int(num)))
                
                t0 = time.perf_counter()
                if key in map_u:
                    mea = map_u[key]
                    old_value = mea.meas_value_ideal
                    run_dpdse.update_measurement(str(uuid), measurement.MeasType(int(num)), value, map_u, value_in_pu=False, data='simulation') 
                    if uuid == '_4C65164895AF4C7F8C748A4A65F70CD5':
                         if int(num) == 7:
                            print(f"---{i}---")
                            v_tillmannshof_mag[i] = mea.meas_value_act # measured voltage mag at fault bus
                         elif int(num) == 8:
                            v_tillmannshof_ang[i] = mea.meas_value_act # measured voltage ang at fault bus
                elif key in map_z:
                    mea = map_z[key]
                    old_value = mea.meas_value_ideal
                    run_dpdse.update_measurement(str(uuid), measurement.MeasType(int(num)), value, map_z, value_in_pu=False, data='simulation') 
                else:
                    pass
            except ValueError:
                pass  # skip malformed column names
    
    run_dpdse.predict() 
    run_dpdse.correct()

    run_dpdse.prepare_output(out='est') # output based on estimates

    # Perform state estimation
    #state_estimation_results = nv_state_estimator.DsseCall(system, measurements_set)

    
    # Extract values from the dictionary and append to the row
    for key, val in run_dpdse.states_output.items():
        #print("output: ", key, val)
        for array in val:
            each_row.append(array.item()) # Get values (real, imag, mag, phase)
    est_values.append(each_row)
    #print("each row: ", np.shape(each_row))
    
    # convert noisy measurement of mag and angle to complex
    v_tillmannshof_re = v_tillmannshof_mag * np.cos(v_tillmannshof_ang)
    v_tillmannshof_im = v_tillmannshof_mag * np.sin(v_tillmannshof_ang)

    # extract WLS results
    #for br in state_estimation_results.branches:
    #    if br.topology_branch.uuid == '_9764573B76AF49F3AB8434F2F9A3897F':
    #       i_ltg_l07H_t02h_re_wls[i] = np.real(br.current)
    #        i_ltg_l07H_t02h_im_wls[i] = np.imag(br.current)
    i = i + 1


output_df = pd.DataFrame(est_values, columns=columns_names)
output_df["v_tillmannshof_re_meas"] = v_tillmannshof_re
output_df["v_tillmannshof_im_meas"] = v_tillmannshof_im
output_df['i_ltg_l07H_t02h_im_wls'] = i_ltg_l07H_t02h_im_wls
output_df['i_ltg_l07H_t02h_re_wls'] = i_ltg_l07H_t02h_re_wls
output_df.to_csv('Seguro_net2_DSEoutput_ts_large.csv', index=False)

# Put file into storage
#store.put_file("data/results/state-estimation/dynamic/Seguro_net2_DSEoutput_ts_large.csv", "Seguro_net2_DSEoutput_ts_large.csv")


print("DP-DSE FOR SEGURO NET2 COMPLETE!!")

