import numpy as np
from enum import Enum
from pyvolt import measurement
from pyvolt import network as net
from pyvolt import nv_powerflow
from pyvolt import nv_state_estimator
import scipy as spy
import observability as ob
import pandas as pd
from scipy.signal import StateSpace, cont2discrete


class Line_Type(Enum):
    RL = 1
    PI = 2


def phasor_complex(mag, ph):
    # mag and ph are list of magnitudes and phase angles
    # returns an array of complex values
    if len(mag) == 0 or len(ph) == 0:
        raise Exception("Empty list. Nothing to convert")
    elif len(mag) != len(ph):
        raise Exception("Magnitude and phase list lengths are not same. There is either a missing magnitude or a "
                        "phase angle")
    else:
        m_arr = np.array(mag)
        ph_arr = np.array(ph)
        cmplx = m_arr * np.exp(1j * ph_arr)
        return cmplx


class DpDse:
    '''
    DpDse Class 
    SV: all load node voltages and all branch currents
    control inputs: all load injection (currents or power for PI model, only power for RL model) and all generator voltages
    '''
    # create A and B matrix for network ,
    # PI: SV: branch currents and load voltages in rectangular form
    # PI: Both banch currents and load voltages are differential equations
    # PI: u : generator voltages and load current injections in rectangular form (curent injections can be obtained from direct measuerments of injection or load powers)
    # either load current measurements are available or power injection measurements are available
    # RL: branch currents and load voltages and load current injections in rectangular form
    # RL: Only branch currents have differential equations, load voltages are calculated as algebraic equations
    # RL: u : generator voltages and power injections 
    def __init__(self, network, measurement_set, time_step, line_type=Line_Type.RL):
        if not isinstance(network, net.System):
            raise Exception("network must be an object of class Network of PyVolt")

        if not isinstance(measurement_set, measurement.MeasurementSet):
            raise Exception("measurement_set must be an object of class MeasurementSet of PyVolt")

        if not isinstance(line_type, Line_Type):
            raise Exception("line_type must be an object of class Line_Type")

        self.network = network
        self.measurement_set = measurement_set
        self.time_step = time_step
        self.line_type = line_type
        self.Adt = np.empty((0, 0))  # descretized A matrix
        self.Bdt = np.empty((0, 0))  # descretized B matrix
        self.Act = np.empty((0, 0))  # continuous A matrix
        self.Bct = np.empty((0, 0))  # continuous B matrix
        self.num_sv = 0  # number of state variables
        self.num_u = 0  # number of control inputs
        self.num_z = 0  # number of measurements (not including control inputs)
        self.num_g = 0  # number of generator nodes
        self.num_l = 0  # number of load nodes; including zero injection nodes
        self.num_b = 0  # number of branches
        self.u = measurement.MeasurementSet()  # list of control inputs in the form [vmag, vph, P, Q] or [vmag, vp, Iinj_mag, Iinj_ph] (TODO: second one needs new measurement type current injection in PyVolt)
        self.x_est = np.array([])  # list of values of estimated states
        self.x_pred = np.array([])  # list of values of predicted states
        self.P_est = np.array([])  # estimation covariance
        self.P_pred = np.array([])  # prediction covariance
        self.z = measurement.MeasurementSet()  # array of measurement objects which will be used for correction step
        self.R = np.empty((0, 0)) # measurement covariance matrix
        self.V = np.empty((0, 0)) # control input covariance matrix
        self.states_output = dict()
        self.tc = 100 # time constant for load voltage calculations

    def initialize_dse(self):
        fo = 50  # TODO: where can we get network nominal frequency information
        w_o = 2 * np.pi * fo

        # set type of line and if load resistance is to be considered # TODO: How and where to specify load resistance?
        load_resistance = False

        # TODO: How to set P_rLoad and R_L??
        # Set base quantities
        Vbase = 12.66  # line-line voltage #TODO: How to set vbase?
        base_apparent_power = 1  # MVA

        get_ES_node_index = [gen_node.index for gen_node in self.network.get_ES_nodes()]
        get_EC_node_index = [load_node.index for load_node in self.network.get_EC_nodes()]

        num_gen = len(get_ES_node_index)
        num_load = len(get_EC_node_index)
        num_nodes = num_gen + num_load
        num_branch = self.network.get_branch_num()

        self.num_b = num_branch
        self.num_g = num_gen
        self.num_l = num_load

        P_rLoad = np.array(
            [[1], [.9], [0.6], [0.6], [0.6], [1.5]])  # this can be built from the power flow results or CIM SV file
        P_rLoad = np.ones((num_load, 1)) # TODO: not really correct!

        R_L = np.divide(Vbase * Vbase * np.ones((num_load, 1)), P_rLoad).flatten()
        R_L = np.diag(R_L)

        # Form bus-branch incidence matrix
        PS_A = self.network.get_bus_branch_incidence_matrix()
        # Separating A_G and A_L
        PS_A_G = PS_A[get_ES_node_index, :]
        PS_A_L = PS_A[get_EC_node_index, :]
        self.Al = PS_A_L

        # Extracting R, L, C
        branch_r = [1e-15 if x == 0.0 else x for x in self.network.get_branch_R()]
        cables_R = np.diag(branch_r)
        branch_x = [1e-15 if x == 0.0 else x for x in
                    self.network.get_branch_X()]
        cables_L = np.diag(branch_x) / w_o
        cables_C = np.diag(np.dot(abs(PS_A_L), self.network.get_branch_BCH()) * 0.5) / w_o

        # create A, B matrices
        # creating A matrix
        A11 = np.dot(-np.linalg.inv(cables_L), cables_R)
        A12 = w_o * np.eye(num_branch)
        A13 = np.dot(np.linalg.inv(cables_L), PS_A_L.transpose())
        A14 = np.zeros((num_branch, num_load))

        A21 = -w_o * np.eye(num_branch)
        A22 = np.dot(-np.linalg.inv(cables_L), cables_R)
        A23 = np.zeros((num_branch, num_load))
        A24 = np.dot(np.linalg.inv(cables_L), PS_A_L.transpose())

        if self.line_type == Line_Type.PI:
            A31 = np.dot(-np.linalg.inv(cables_C), PS_A_L)
        else:
            A31 = np.zeros((num_load, num_branch))
        A32 = np.zeros((num_load, num_branch))
        if self.line_type == Line_Type.PI and load_resistance is True:
            A33 = -np.linalg.inv(np.dot(cables_C, R_L))
        elif self.line_type == Line_Type.PI and load_resistance is False:
            A33 = np.zeros((num_load, num_load))  # when load resistance is not to be considered TODO: Is this correct? Even for PI model when load resis is not there? 
        elif self.line_type == Line_Type.RL: # TODO: is this correct for both with and without load resis? If load resistance is there, shouldnt it change? 
            A33 = -(1/self.tc)*np.eye(num_load) 
        if self.line_type == Line_Type.PI:
            A34 = w_o * np.eye(num_load)
        else:
            A34 = np.zeros((num_load, num_load))  # when shunt capacitance is not to be considered

        A41 = np.zeros((num_load, num_branch))
        if self.line_type == Line_Type.PI:
            A42 = np.dot(-np.linalg.inv(cables_C), PS_A_L)
        else:
            A42 = np.zeros((num_load, num_branch))
        if self.line_type == Line_Type.PI:
            A43 = -w_o * np.eye(num_load)
        else:
            A43 = np.zeros((num_load, num_load))  # when shunt capacitance is not to be considered
        if self.line_type == Line_Type.PI and load_resistance is True:
            A44 = -np.linalg.inv(np.dot(cables_C, R_L))
        elif self.line_type == Line_Type.PI and load_resistance is False:
            A44 = np.zeros((num_load, num_load))  # when load resistance is not to be considered
        elif self.line_type == Line_Type.RL:
            A44 = -(1/self.tc)*np.eye(num_load)

        # Creating B matrix
        B11 = np.dot(np.linalg.inv(cables_L), PS_A_G.transpose())
        B12 = np.zeros((num_branch, num_gen))
        B13 = np.zeros((num_branch, num_load))
        B14 = np.zeros((num_branch, num_load))

        B21 = np.zeros((num_branch, num_gen))
        B22 = np.dot(np.linalg.inv(cables_L), PS_A_G.transpose())
        B23 = np.zeros((num_branch, num_load))
        B24 = np.zeros((num_branch, num_load))

        B31 = np.zeros((num_load, num_gen))
        B32 = np.zeros((num_load, num_gen))
        if self.line_type == Line_Type.PI:
            B33 = -np.linalg.inv(cables_C)
        else:
            B33 = (1/self.tc)*np.eye(num_load)
        B34 = np.zeros((num_load, num_load))

        B41 = np.zeros((num_load, num_gen))
        B42 = np.zeros((num_load, num_gen))
        B43 = np.zeros((num_load, num_load))
        if self.line_type == Line_Type.PI:
            B44 = -np.linalg.inv(cables_C)
        else:
            B44 = (1/self.tc)*np.eye(num_load)


        self.num_sv = 2 * num_branch + 2 * num_load
        self.num_u = 2 * num_gen + 2 * num_load

        SS_A = np.bmat([[A11, A12, A13, A14], [A21, A22, A23, A24], [A31, A32, A33, A34], [A41, A42, A43, A44]])
        SS_B = np.bmat([[B11, B12, B13, B14], [B21, B22, B23, B24], [B31, B32, B33, B34], [B41, B42, B43, B44]])
        self.Act = SS_A
        self.Bct = SS_B
        
        Ad, Bd, Cd, Dd, _ = cont2discrete(( SS_A, SS_B, np.eye(self.num_sv), np.zeros((self.num_sv, self.num_u)) ), self.time_step, method = 'zoh')
        self.Adt = Ad
        self.Bdt = Bd

        # map sv index to uuids
        self.set_sv_idx_uuid()

        # prepare dictionary to store output
        self.set_states_output_dict()

        # initialize initial estimation covariance
        self.P_est = 1e-10 * np.eye(self.num_sv)

        # initialize process error
        self.Q = 1e-10 * np.eye(self.num_sv)

       
        # initialize SV using static se or power flow
        try:
            # TODO: the current injection as measurement is not considered in PyVolt! Pyvolt will throw error!
            static_se_results = nv_state_estimator.DsseCall(self.network, self.measurement_set)
            self.initialize_sv(static_se_results)
        except:
            print("PyVolt SE threw error! State variables are initialized from power flow results!")
            static_se_results, num_iter = nv_powerflow.solve(self.network)
            self.initialize_sv(static_se_results)
        
        # extract the control variables and measurements separately
        self.separate_inputs()

        # initialize measurement error
        self.R = 1e-10 * np.eye(self.num_z)

        # initialize control input measurement error
        self.V = 1e-10 * np.eye(self.num_u)

        self.assemble_u()

        # check observability and detectability of available measurement configuration
        _, H = self.build_hx_H() 
        is_ob, O = ob.is_system_observable(self.Adt, H)
        is_detect = ob.is_system_detectable(self.Adt, H, system_type='Discrete')
        v = ob.recommended_sv_measurements(self.Adt, H)
        print("Is system observable: ", is_ob)
        print("Is system detectable: ", is_detect)

    def separate_inputs(self):
        # from the entire measurement_set, extract and separate the measurements which form control variables
        # (vg, il or sl) and rest of the measurements which form observations used for correct step

        gen_uuid = [gen_node.uuid for gen_node in self.network.get_ES_nodes()]
        load_uuid = [load_node.uuid for load_node in self.network.get_EC_nodes()]

        # first extract all measurements of a node of type voltage and power injections of generator nodes and load
        # nodes respectively, the remaining goes into z list
        gen_input_type = [measurement.MeasType.Vpmu_mag, measurement.MeasType.Vpmu_phase]
        load_power_input_type = [measurement.MeasType.Sinj_real, measurement.MeasType.Sinj_imag]
        load_current_input_type = [measurement.MeasType.Ipmu_inj_mag, measurement.MeasType.Ipmu_inj_phase]

        for meas in self.measurement_set.measurements:
            if meas.element_type == measurement.ElemType.Node:
                if meas.meas_type in gen_input_type and meas.element.uuid in gen_uuid:
                    self.u.measurements.append(meas)
                if meas.meas_type in load_power_input_type and meas.element.uuid in load_uuid:
                    meas.meas_value_act = meas.meas_value_act 
                    meas.meas_value = meas.meas_value
                    self.u.measurements.append(meas)
                if meas.meas_type in load_current_input_type and meas.element.uuid in load_uuid and self.line_type == Line_Type.PI: # only PI model can take current injection as input, in RL model, current injection is measurement
                    if meas.meas_type == measurement.MeasType.Ipmu_inj_phase:
                        meas.meas_value_act = meas.meas_value_act  
                        meas.meas_value = meas.meas_value 
                    self.u.measurements.append(meas)
                    
        # sort control inputs in required order 
        sorted_u = measurement.MeasurementSet()

        # Sort measurements  in the order required by the SE algorithm
        # Required order: Vpmu_mag, Vpmu_phase - in the respective order of gen, [Sinj_real, Sinj_imag or Ipmu_inj_mag, Ipmu_inj_phase] - in the respective order of load
        u_uid = gen_uuid + load_uuid
        for uid in u_uid:
            for meas in self.u.measurements:
                if meas.element.uuid == uid: 
                    sorted_u.measurements.append(meas)
        
        self.u = sorted_u
    
        # extract z
        for item in self.measurement_set.measurements:
            if item not in self.u.measurements:
               self.z.measurements.append(item)
        
        # sort z into required order
        self.z = self.z.getSortedMeasurementSet() # sort measurements by type

        # set number of measurements and control inputs
        self.num_z = len(self.z.measurements)
        self.num_u = len(self.u.measurements)
        

    def update_covariance_pmu(self, cov, index_mag, index_phase, u_z_type="z"):
        # find covariance of phasor measurements in rectangular
        if len(index_mag) != 0 and len(index_phase) != 0:
            for index, (idx_mag, idx_theta) in enumerate(zip(index_mag, index_phase)):
                if u_z_type == 'z':
                    value_amp = self.z.measurements[idx_mag].meas_value_act
                    value_theta = self.z.measurements[idx_theta].meas_value_act
                elif u_z_type == 'u':
                    value_amp = self.u.measurements[idx_mag].meas_value_act # TODO: this may be wrong! Because, in assemble_u the getActuals is called, which already replaces mag and phase act values with real and imag! 
                    value_theta = self.u.measurements[idx_theta].meas_value_act
                rot_mat = np.array([[np.cos(value_theta), - value_amp * np.sin(value_theta)],
                                    [np.sin(value_theta), value_amp * np.cos(value_theta)]])
                starting_cov = np.array([[cov[idx_mag], 0], [0, cov[idx_theta]]])
                final_cov = np.inner(rot_mat, np.inner(starting_cov, rot_mat.transpose()))
                if u_z_type == 'z':
                    self.R[idx_mag][idx_mag] = final_cov[0][0]
                    self.R[idx_theta][idx_theta] = final_cov[1][1]
                    self.R[idx_mag][idx_theta] = final_cov[0][1]
                    self.R[idx_theta][idx_mag] = final_cov[1][0]
                elif u_z_type == 'u':
                    self.V[idx_mag][idx_mag] = final_cov[0][0]
                    self.V[idx_theta][idx_theta] = final_cov[1][1]
                    self.V[idx_mag][idx_theta] = final_cov[0][1]
                    self.V[idx_theta][idx_mag] = final_cov[1][0]

    
    def assemble_u(self):
        u = []
        # the current injection from load power to be calculated only for PI model. For RL, direct power can be used
        sinj_real_idx = self.u.getIndexOfMeasurements(measurement.MeasType.Sinj_real)
        sinj_imag_idx = self.u.getIndexOfMeasurements(measurement.MeasType.Sinj_imag)
        ib_est_re = self.x_est[0 : self.num_b]
        ib_est_im = self.x_est[self.num_b : 2*self.num_b]
        il_inj_re = self.Al @ ib_est_re
        il_inj_im = self.Al @ ib_est_im
        
        # also prepare the control input covariance
        u_covar = self.u.getCovarianceMatrixActuals()
        self.V = np.diag(u_covar)

        # uncertainty propagation from phasor to complex measurements of PMU - updating measurement covariance
        u_Il_mag_idx = self.u.getIndexOfMeasurements(measurement.MeasType.Ipmu_inj_mag)
        u_Il_phase_idx = self.u.getIndexOfMeasurements(measurement.MeasType.Ipmu_inj_phase)
        u_Vl_mag_idx = self.u.getIndexOfMeasurements(measurement.MeasType.Vpmu_mag)
        u_Vl_phase_idx = self.u.getIndexOfMeasurements(measurement.MeasType.Vpmu_phase)
        self.update_covariance_pmu(u_covar, u_Il_mag_idx, u_Il_phase_idx, u_z_type='u')
        self.update_covariance_pmu(u_covar, u_Vl_mag_idx, u_Vl_phase_idx, u_z_type='u')

        for index, (idx_re, idx_im) in enumerate(zip(sinj_real_idx, sinj_imag_idx)):
            p = (-1) * self.u.measurements[idx_re].meas_value_act # negating to make it injection convention
            q = (-1) * self.u.measurements[idx_im].meas_value_act
            if self.u.measurements[idx_re].element.uuid != self.u.measurements[idx_im].element.uuid:
                print("Real and reactive power do not belong to same element!")
            else:
                uid = self.u.measurements[idx_re].element.uuid
                v_est_re = self.x_est[self.vl_re_idx_uuid[uid], 0]
                v_est_im = self.x_est[self.vl_im_idx_uuid[uid], 0]
                if self.line_type == Line_Type.PI:
                    v_sq = (v_est_re * v_est_re + v_est_im * v_est_im)
                    il_re = (v_est_re * p + v_est_im * q) / v_sq
                    il_im = (v_est_im * p - v_est_re * q) / v_sq
                    #print(f"uuid: {uid}, p: {p}, q: {q}, il_re: {il_re}, il_im: {il_im}")
                    self.u.measurements[idx_re].meas_value_act = il_re
                    self.u.measurements[idx_im].meas_value_act = il_im
                    # compute covariance uncertainty propagation for Pl Ql --> Il
                    rot_mat = np.array([[v_est_re/v_sq , v_est_im/v_sq],
                                        [v_est_im/v_sq, - v_est_re/v_sq]])
                    starting_cov = np.array([[u_covar[idx_re], 0], [0, u_covar[idx_im]]])
                    final_cov = np.inner(rot_mat, np.inner(starting_cov, rot_mat.transpose()))
                    self.V[idx_re][idx_re] = final_cov[0][0]
                    self.V[idx_im][idx_im] = final_cov[1][1]
                    self.V[idx_re][idx_im] = final_cov[0][1]
                    self.V[idx_im][idx_re] = final_cov[1][0]

                elif self.line_type == Line_Type.RL:
                    load_node_index = self.vl_re_idx_uuid[uid] - 2*self.num_b # the imag index is same as real index for the use here, because il is divided into real and imag il_inj
                    if np.isclose(il_inj_re[load_node_index], 0, atol=1e-10) and np.isclose(il_inj_im[load_node_index], 0, atol=1e-10): 
                        self.u.measurements[idx_re].meas_value_act = -v_est_re # negate to ensure the injection convention
                        self.u.measurements[idx_im].meas_value_act = -v_est_im
                    else:
                        il_sq = (il_inj_re[load_node_index] * il_inj_re[load_node_index] + il_inj_im[load_node_index] * il_inj_im[load_node_index])
                        self.u.measurements[idx_re].meas_value_act = (p * il_inj_re[load_node_index] - q * il_inj_im[load_node_index])/il_sq
                        self.u.measurements[idx_im].meas_value_act = (q * il_inj_re[load_node_index] + p * il_inj_im[load_node_index])/il_sq
                        
                        # compute covariance uncertainty propagation for Pl Ql --> Vl
                        rot_mat = np.array([[(il_inj_re[load_node_index]/il_sq).item() , - (il_inj_im[load_node_index]/il_sq).item()],
                                            [(il_inj_im[load_node_index]/il_sq).item() , (il_inj_re[load_node_index]/il_sq).item()]])
                        starting_cov = np.array([[u_covar[idx_re], 0], [0, u_covar[idx_im]]])
                        final_cov = np.inner(rot_mat, np.inner(starting_cov, rot_mat.transpose()))
                        self.V[idx_re][idx_re] = final_cov[0][0]
                        self.V[idx_im][idx_im] = final_cov[1][1]
                        self.V[idx_re][idx_im] = final_cov[0][1]
                        self.V[idx_im][idx_re] = final_cov[1][0]                      

        # the current injections (in PI) and the generator voltages is now converted to real-imag 
        u = self.u.getMeasValuesActuals().reshape((-1, 1))

        all_u_index = []
        all_u_index.extend(u_Vl_mag_idx)
        all_u_index.extend(u_Vl_phase_idx)
        all_u_index.extend(sinj_real_idx)
        all_u_index.extend(sinj_imag_idx)
        #u_new = np.array([self.u.measurements[m].meas_value_act.item() if isinstance(self.u.measurements[m].meas_value_act, np.ndarray) else self.u.measurements[m].meas_value_act   for m in all_u_index] ).reshape(-1, 1)
        u_new = np.array([u[m] for m in all_u_index] ).reshape(-1, 1)
        #print("assemble_u u: ", u, np.shape(u)) # TODO: the above only separates Vl and S, but doesnt account Il (this is required for PI) - yet to do!
        #print("assemble_u u_new: ", u_new, np.shape(u_new))

        return u_new
        
   

    def initialize_sv(self, static_se_results):
        # Initialize state variables (in actuals, per phase)
        ib_re = []
        ib_im = []
        vl_re = []
        vl_im = []
        for br in static_se_results.branches:
            ib_re.append(br.current.real)
            ib_im.append(br.current.imag)
        for node in static_se_results.nodes:
            if node.topology_node.type == net.BusType.PQ:
                vl_re.append(node.voltage.real)
                vl_im.append(node.voltage.imag)
        x_init = np.concatenate((ib_re, ib_im, vl_re, vl_im)).reshape((-1, 1))
        self.x_est = x_init 
        self.x_pred = x_init # x_pred is initialized to enable observability check before the estimator begins.
      

    def predict(self):
        # assemble u vector into required form and compute covariance
        u = self.assemble_u()     

        # predict the states for next time-step
        self.x_pred = self.Adt @ self.x_est + self.Bdt @ u

        # compute prediction covariance
        self.P_pred = self.Adt @ self.P_est @ (self.Adt).T + self.Bdt @ self.V @ (self.Bdt).T + self.Q

    def correct(self):
        # Step 1: Build measurement covariance matrix
        z_covar = self.z.getCovarianceMatrixActuals()
        self.R = np.diag(z_covar)
        
        # Step 2: uncertainty propagation from phasor to complex measurements of PMU - updating measurement covariance
        z_Ib_mag_idx = self.z.getIndexOfMeasurements(measurement.MeasType.Ipmu_mag)
        z_Ib_phase_idx = self.z.getIndexOfMeasurements(measurement.MeasType.Ipmu_phase)
        z_Vl_mag_idx = self.z.getIndexOfMeasurements(measurement.MeasType.Vpmu_mag)
        z_Vl_phase_idx = self.z.getIndexOfMeasurements(measurement.MeasType.Vpmu_phase)
        self.update_covariance_pmu(z_covar, z_Ib_mag_idx, z_Ib_phase_idx, u_z_type='z')
        self.update_covariance_pmu(z_covar, z_Vl_mag_idx, z_Vl_phase_idx, u_z_type='z')
        
        # Step 3: get all measurements, phasor values are subsituted by complex
        z = np.array(self.z.getMeasValuesActuals()).reshape((-1,1))

       
        # Step 4: Build measurement functions (h_x) and Jacobians (H)        
        h_x, H = self.build_hx_H()       
        
        # Step 5: compute Kalman gain
        S = np.linalg.inv(H @ self.P_pred @ H.T + self.R) 
        K = self.P_pred @ H.T @ S

        # Step 6: calculate state estimates
        self.x_est = self.x_pred + K @ (z - h_x)

        # Step 7. compute estimation covariance
        self.P_est = (np.eye(self.num_sv) - K @ H) @ self.P_pred
               
        # preparing and updating dictionary structure for output: {(uid, type) : (real, imag, mag, phase, real_var, imag_var)} 
        # Convert dictionary to pandas Series
        state_series = pd.Series(self.states_output)
        # Update all values by directly assigning the updated list
        state_series[:] = self.prepare_output()
        # Convert Series back to dictionary
        self.states_output = state_series.to_dict()


        return 1

    def build_hx_H(self):   
        # Required order of measurements: Vmag, Sinj_real, Sinj_imag, S1_real, S1_imag, S2_real, S2_imag, Imag, Vpmu_mag, Vpmu_phase, Ipmu_mag, Ipmu_phase, Ipmu_inj_mag, Ipmu_inj_phase
        # [Sinj_real, Sinj_imag] - are part of control inputs (u), and dont appear as measurements (z)  
        
        # build h and H for load voltage phasors [Vpmu_mag, Vpmu_phase]
        z_Vl_mag_idx = self.z.getIndexOfMeasurements(measurement.MeasType.Vpmu_mag)
        z_Vl_phase_idx = self.z.getIndexOfMeasurements(measurement.MeasType.Vpmu_phase)
        h1, H1, h2, H2 = self.hx_H_load_voltages_phasor(len(z_Vl_mag_idx), z_Vl_mag_idx, z_Vl_phase_idx)
        
        # build h and H for branch current phasors [Ipmu_mag, Ipmu_phase]
        z_Ib_mag_idx = self.z.getIndexOfMeasurements(measurement.MeasType.Ipmu_mag)
        z_Ib_phase_idx = self.z.getIndexOfMeasurements(measurement.MeasType.Ipmu_phase)
        h3, H3, h4, H4 = self.hx_H_branch_current_phasor(len(z_Ib_mag_idx), z_Ib_mag_idx, z_Ib_phase_idx)

        # build h and H for branch current magnitude [I_mag]
        ib_idx_mag = self.z.getIndexOfMeasurements(measurement.MeasType.I_mag)
        h5, H5 = self.hx_H_branch_current_magnitude(len(ib_idx_mag), ib_idx_mag)
        
        # build h and H for branch powers [S1_real, S1_imag, S2_real, S2_imag]

        # build h and H for load voltage magnitudes [V_mag]
        vl_idx_mag = self.z.getIndexOfMeasurements(measurement.MeasType.V_mag)
        h6, H6 = self.hx_H_load_voltage_magnitude(len(vl_idx_mag), vl_idx_mag)

        # build h and H for load and generator current injections [Ipmu_inj_mag, Ipmu_inj_phase]
        
        # stack Jacobians 
        H = np.concatenate((H6, H5, H1, H2, H3, H4), axis=0)
       
        # stack measurement functions
        h_x = np.concatenate((h6, h5, h1, h2, h3, h4), axis=0) 

        return h_x, H
    
    
    def hx_H_load_voltages_phasor(self, nvl, index_vre, index_vim):
        # at every iteration we update h(x) vector where Vl measure are available
        h1 = np.zeros((nvl, 1))
        h2 = np.zeros((nvl, 1))
        # the Jacobian rows where voltage measurements are presents is updated
        H1 = np.zeros((nvl, self.num_sv))
        H2 = np.zeros((nvl, self.num_sv))

        for i, (idx_vre, idx_vim) in enumerate(zip(index_vre, index_vim)):
            # get index of the node
            node_uuid_re = self.z.measurements[idx_vre].element.uuid
            node_uuid_im = self.z.measurements[idx_vim].element.uuid
            if node_uuid_re == node_uuid_im:
                sv_re_idx = self.vl_re_idx_uuid[node_uuid_re]
                sv_im_idx = self.vl_im_idx_uuid[node_uuid_re]
                h1[i][0] = self.x_pred[sv_re_idx]
                H1[i][sv_re_idx] = 1
                h2[i][0] = self.x_pred[sv_im_idx]
                H2[i][sv_im_idx] = 1
            else:
                print("Real and imaginary for load voltage do not belong to same element!")
        return h1, H1, h2, H2
    
    def hx_H_branch_current_phasor(self, nib, index_ibre, index_ibim):
        # at every iteration we update h(x) vector where Ib measure are available
        h3 = np.zeros((nib, 1))
        h4 = np.zeros((nib, 1))
        # the Jacobian rows where branch current measurements are presents is updated
        H3 = np.zeros((nib, self.num_sv))
        H4 = np.zeros((nib, self.num_sv))

        for i, (idx_ibre, idx_ibim) in enumerate(zip(index_ibre, index_ibim)):
            # get index of the node
            node_uuid_re = self.z.measurements[idx_ibre].element.uuid
            node_uuid_im = self.z.measurements[idx_ibim].element.uuid
            if node_uuid_re == node_uuid_im:
                sv_re_idx = self.ib_re_idx_uuid[node_uuid_re]
                sv_im_idx = self.ib_im_idx_uuid[node_uuid_re]
                h3[i][0] = self.x_pred[sv_re_idx]
                H3[i][sv_re_idx] = 1
                h4[i][0] = self.x_pred[sv_im_idx]
                H4[i][sv_im_idx] = 1
            else:
                print("Real and imaginary for branch currents do not belong to same element!")
        return h3, H3, h4, H4
    
    def hx_H_branch_current_magnitude(self, nib, index_ib_mag):
        # at every iteration we update h(x) vector where Ib measure are available
        h5 = np.zeros((nib, 1))

        # the Jacobian rows where branch current measurements are presents is updated
        H5 = np.zeros((nib, self.num_sv))

        for i, idx_ib in enumerate(index_ib_mag):
            # get index of the node
            node_uuid = self.z.measurements[idx_ib].element.uuid
            sv_re_idx = self.ib_re_idx_uuid[node_uuid]
            sv_im_idx = self.ib_im_idx_uuid[node_uuid]
            h5[i][0] = np.sqrt(self.x_pred[sv_re_idx] ** 2 + self.x_pred[sv_im_idx] ** 2)
            H5[i][sv_re_idx] = self.x_pred[sv_re_idx] / np.sqrt(self.x_pred[sv_re_idx] ** 2 + self.x_pred[sv_im_idx] ** 2)
            H5[i][sv_im_idx] = self.x_pred[sv_im_idx] / np.sqrt(self.x_pred[sv_re_idx] ** 2 + self.x_pred[sv_im_idx] ** 2)

        return h5, H5
    
    def hx_H_load_voltage_magnitude(self, nvl, index_vl_mag):
        # at every iteration we update h(x) vector where Ib measure are available
        h6 = np.zeros((nvl, 1))

        # the Jacobian rows where branch current measurements are presents is updated
        H6 = np.zeros((nvl, self.num_sv))

        for i, idx_vl in enumerate(index_vl_mag):
            # get index of the node
            node_uuid = self.z.measurements[idx_vl].element.uuid
            sv_re_idx = self.vl_re_idx_uuid[node_uuid]
            sv_im_idx = self.vl_im_idx_uuid[node_uuid]
            h6[i][0] = np.sqrt(self.x_pred[sv_re_idx] ** 2 + self.x_pred[sv_im_idx] ** 2)
            H6[i][sv_re_idx] = self.x_pred[sv_re_idx] / np.sqrt(self.x_pred[sv_re_idx] ** 2 + self.x_pred[sv_im_idx] ** 2)
            H6[i][sv_im_idx] = self.x_pred[sv_im_idx] / np.sqrt(self.x_pred[sv_re_idx] ** 2 + self.x_pred[sv_im_idx] ** 2)

        return h6, H6


    def update_measurement(self, element_uuid, meas_type, meas_data, mapping, value_in_pu=True):
        """
        to update meas_value of a specific measurement object in the measurements array
        """
        # TODO: How to know if node current and node power has positive value for consumption (out-of grid) or injection (into-grid)???
        # In PyVolt program, it is assumed that node current and node power has positive if it is consumption (out-of grid), therefore negated to make it injection - in assemble_u function of dpdse.py
        # From RTDS, we get directly the injection currents, therefore not negated here. 
        #TODO: Need a standard way to fix this.
        
        load_power_input_type = [measurement.MeasType.Sinj_real, measurement.MeasType.Sinj_imag]
        load_current_input_type = [measurement.MeasType.Ipmu_inj_mag, measurement.MeasType.Ipmu_inj_phase]

        key = (element_uuid, meas_type)
        if key in mapping:
            meas = mapping[key]
        # only update measurements that are already included in the measurements set
            if meas.element.uuid == element_uuid:
                # special case for voltage magnitude: SOGNO interface only knows Vpmu_mag while measurement set distincts between Vpmu_mag and V_mag
                if (meas.meas_type == measurement.MeasType.Vpmu_mag or meas.meas_type == measurement.MeasType.V_mag):
                    # volt pu conversion assuming that meas_data from device are in volts and single-phase value according to sogno interface
                    # while baseVoltage from CIM in [kV] and three-phase value
                    if not value_in_pu:
                        meas_value_pu = meas_data / (meas.element.baseVoltage / np.sqrt(3) )
                        meas_value_act = meas_data
                    else:
                        meas_value_pu = meas_data
                        meas_value_act = meas_data * (meas.element.baseVoltage / np.sqrt(3) )
                    #print("Updating measurement value for {} of type {} from {} to {}, value_in_pu {} ".format(meas.element.uuid, str(meas.meas_type), meas.meas_value_act, meas_value_act, value_in_pu))                    
                    meas.meas_value = meas_value_pu 
                    meas.meas_value_act = meas_value_act                
                # special case for voltage magnitude: SOGNO interface only knows Ipmu_mag while measurement set distincts between Ipmu_mag and I_mag
                elif (meas.meas_type == measurement.MeasType.Ipmu_mag or meas.meas_type == measurement.MeasType.I_mag or meas.meas_type == measurement.MeasType.Ipmu_inj_mag):
                    # current pu conversion assuming that meas_data from device are in A and single-phase value according to sogno interface
                    # while base_current in [kA]
                    if not value_in_pu:
                        meas_value_pu = meas_data / (meas.element.base_current )
                        meas_value_act = meas_data
                    else:
                        meas_value_pu = meas_data
                        meas_value_act = meas_data * (meas.element.base_current )
                    #print("Updating measurement value for {} of type {} from {} to {}, value_in_pu {}".format(meas.element.uuid, str(meas.meas_type), meas.meas_value_act, meas_value_act, value_in_pu))
                    meas.meas_value = meas_value_pu
                    meas.meas_value_act = meas_value_act
                # case for other measurements 
                elif (meas_type == meas.meas_type and (meas_type == measurement.MeasType.S1_real or meas_type == measurement.MeasType.S1_imag or meas_type == measurement.MeasType.Sinj_real or meas_type == measurement.MeasType.Sinj_imag)): 
                    # power pu conversion assuming that meas_data from device are in watts and single-phase value according to sogno interface
                    # while baseApparent power in [MW] and three-phase value
                    if not value_in_pu:    
                        if (meas_type == measurement.MeasType.Sinj_real or meas_type == measurement.MeasType.Sinj_imag):
                            meas_value_pu = meas_data / (meas.element.base_apparent_power / 3 )
                            meas_value_act = meas_data
                        else:
                            meas_value_pu = meas_data / (meas.element.base_apparent_power / 3 )
                            meas_value_act = meas_data
                    else:
                        if (meas_type == measurement.MeasType.Sinj_real or meas_type == measurement.MeasType.Sinj_imag): # injection power should be negated if consumption is positive 
                            meas_value_pu = meas_data
                            meas_value_act = meas_data * (meas.element.base_apparent_power / 3 )
                        else:
                            meas_value_pu = meas_data
                            meas_value_act = meas_data * (meas.element.base_apparent_power / 3 )
                    #print("Updating measurement value for {} of type {} from {} to {}, value_in_pu {}".format(meas.element.uuid, str(meas.meas_type), meas.meas_value_act, meas_value_act, value_in_pu))
                    meas.meas_value = meas_value_pu
                    meas.meas_value_act = meas_value_act
                elif (meas_type == meas.meas_type and (meas_type == measurement.MeasType.Vpmu_phase or meas_type == measurement.MeasType.Ipmu_phase or meas_type == measurement.MeasType.Ipmu_inj_phase)):
                    #print("Updating measurement value for {} of type {} from {} to {}, value_in_pu {}".format(meas.element.uuid, str(meas.meas_type), meas.meas_value_act, meas_data, value_in_pu))
                    if meas_type == measurement.MeasType.Ipmu_inj_phase:  
                        meas.meas_value = meas_data    # for negating for current injection
                        meas.meas_value_act = meas_data 
                    else:
                        meas.meas_value = meas_data
                        meas.meas_value_act = meas_data
        else:
            print(f"with uuid={element_uuid} and type={meas_type} not found")




    def set_sv_idx_uuid(self):
        self.ib_re_idx_uuid =  {uuid: index for index, uuid in enumerate(self.getBranchUuid())}
        self.ib_im_idx_uuid =  {uuid: index + self.num_b for index, uuid in enumerate(self.getBranchUuid())}
        self.vl_re_idx_uuid =  {uuid: index + 2*self.num_b for index, uuid in enumerate(self.getLoadUuid())}
        self.vl_im_idx_uuid =  {uuid: index + self.num_l + 2*self.num_b for index, uuid in enumerate(self.getLoadUuid())}
    
    def set_states_output_dict(self):
        curr_dict = {(uuid, 'branch_current'): (0, 0, 0, 0, 0, 0) for index, uuid in enumerate(self.getBranchUuid())}
        volt_dict = {(uuid, 'load_voltage'): (0, 0, 0, 0, 0, 0) for index, uuid in enumerate(self.getLoadUuid())}
        self.states_output = {**curr_dict, **volt_dict}

    def prepare_output(self):
        states = self.x_est 
        covars = np.diag(self.P_est)

        curr_len = self.num_b
        volt_len = self.num_l
        curr_re = np.array(states[:curr_len])
        curr_im = np.array(states[curr_len:2*curr_len])
        volt_re = np.array(states[2*curr_len:2*curr_len + volt_len])
        volt_im = np.array(states[2*curr_len + volt_len:])

        curr_re_var = np.array(covars[:curr_len])
        curr_im_var = np.array(covars[curr_len:2*curr_len])
        volt_re_var = np.array(covars[2*curr_len:2*curr_len + volt_len])
        volt_im_var = np.array(covars[2*curr_len + volt_len:])

        # Compute magnitude
        curr_mag = np.sqrt(curr_re**2 + curr_im**2)
        volt_mag = np.sqrt(volt_re**2 + volt_im**2)
        
        # Compute phase angle (in radians)
        curr_phase = np.arctan2(curr_im, curr_re)
        volt_phase = np.arctan2(volt_im, volt_re)
       
        curr_list = list(zip(curr_re, curr_im, curr_mag, curr_phase, curr_re_var, curr_im_var))
        volt_list = list(zip(volt_re, volt_im, volt_mag, volt_phase, volt_re_var, volt_im_var))
        
        return (curr_list + volt_list)


    def getBranchUuid(self):
        br_uuid = []
        for branch in self.network.branches:
            br_uuid.append(branch.uuid)
        return br_uuid
    
    def getLoadUuid(self):
        l_uuid = []
        for node in self.network.get_EC_nodes():
            l_uuid.append(node.uuid)
        return l_uuid

    def getGenUuid(self):
        g_uuid = []
        for node in self.network.get_ES_nodes():
            g_uuid.append(node.uuid)
        return g_uuid
    
    def extract_branch_currents_estimation(self):
        return self.x_est[:2 * self.num_b]

    def extract_branch_currents_prediction(self):
        return self.x_pred[:2 * self.num_b]

    def extract_load_voltages_estimation(self):
        return self.x_est[2 * self.num_b: 2 * self.num_b + 2 * self.num_l]

    def extract_load_voltages_prediction(self):
        return self.x_pred[2 * self.num_b: 2 * self.num_b + 2 * self.num_l]
    
    def get_num_sv(self):
        return self.num_sv

    def get_num_u(self):
        return self.num_u

    def get_num_g(self):
        return self.num_g

    def get_num_l(self):
        return self.num_l

    def get_num_b(self):
        return self.num_b

    def get_Act(self):
        return self.Act

    def get_Bct(self):
        return self.Bct

    def get_Adt(self):
        return self.Adt

    def get_Bdt(self):
        return self.Bdt
    
    def get_meas_z(self):
        return self.z
    
    def get_meas_u(self):
        return self.u

################################### VALIDATING STATE SPACE CONSTRUCTION #################################
    def check_ss_consistency(self):
        # Execute power flow analysis
        results_pf, num_iter = nv_powerflow.solve(self.network)
        # Print node voltages
        print("Powerflow converged in " + str(num_iter) + " iterations.\n")
        pf_node_voltages = []
        for node in results_pf.nodes:
            pf_node_voltages.append(node.voltage)

        pf_br_currents = []
        for branch in results_pf.branches:
            pf_br_currents.append(branch.current)

        pf_node_currents = [] 
        for node in results_pf.nodes:
            pf_node_currents.append(node.current)

        # initialization of state variables and control inputs
        get_ES_node_index = [gen_node.index for gen_node in self.network.get_ES_nodes()]
        get_EC_node_index = [load_node.index for load_node in self.network.get_EC_nodes()]
        pf_vg = np.array(pf_node_voltages)[np.array(get_ES_node_index)]
        pf_vl = np.array(pf_node_voltages)[np.array(get_EC_node_index)]
        pf_il = np.array(pf_node_currents)[np.array(get_EC_node_index)]
        pf_ibr = np.array(pf_br_currents)
        # TODO: This consistency check is only valid for INPUT_CURRENT_INJECTION model type
        if self.line_type == Line_Type.PI:
            x_init_pf_conv = np.concatenate(
                (pf_ibr.real, pf_ibr.imag, pf_vl.real , pf_vl.imag ))
            u_init_pf_conv = np.concatenate((pf_vg.real , pf_vg.imag, pf_il.real, pf_il.imag))
            # negated load currents as current injections from loads are considered
            u_init_pf_curr_inj = np.concatenate((pf_vg.real , pf_vg.imag , - pf_il.real,
                                                - pf_il.imag))

            print("u_init_pf_curr_inj: ", u_init_pf_curr_inj)
            bu = np.dot(self.Bct, u_init_pf_conv)
            # solve for x = inv(A)*Bu
            bu_curr_inj = np.dot(self.Bct, u_init_pf_curr_inj)
            x_init_from_ss = np.dot(np.linalg.inv(self.Act), -bu_curr_inj.T)

            print("x_init direct pf: ", x_init_pf_conv)
            print("x_init_state space: ", x_init_from_ss)

            print("Ax: ", np.dot(self.Act, x_init_pf_conv))
            print("Bu: ", bu)

            print("Ax_state space: ", np.dot(self.Act, x_init_from_ss).T)
            print("bu_curr_inj: ", bu_curr_inj)
        elif self.line_type == Line_Type.RL:
            x_init_pf_conv = np.concatenate(
                (pf_ibr.real, pf_ibr.imag, pf_vl.real , pf_vl.imag ))
            u_init_pf_conv = np.concatenate((pf_vg.real , pf_vg.imag, pf_vl.real, pf_vl.imag)) # vl_calc is the u here in RL model
            print("u_init_pf_conv : ", u_init_pf_conv, np.shape(u_init_pf_conv) )
            bu = np.dot(self.Bct, u_init_pf_conv)
            
            # inv(A)*B*u will result in singularity!! Therefore only checking if Ax = Bu at steady state.
            #print("Ax: ", np.dot(self.Act, x_init_pf_conv))
            #print("Bu: ", bu)