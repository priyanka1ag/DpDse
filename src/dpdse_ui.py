import numpy as np
from enum import Enum
from pyvolt import measurement
from pyvolt import network as net
from pyvolt import nv_powerflow
from pyvolt import nv_state_estimator
import scipy as spy


# TODO: Only centralized implementation available so far! Feeder selection should be implemented here. Because it affects the state space formation
class DpDse_Centralized(Enum):
    YES = 1
    NO = 2


class DpDse_Model(Enum):
    STANDARD = 1
    AUGMENTED = 2
    INPUT_LOAD_POWER = 3
    LOAD_CURRENT_AUGMENTED = 4


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
    control inputs: all load injection currents and all generator voltages
    '''
    def __init__(self, network, measurement_set, critical_nodes, time_step, model_type=DpDse_Model.STANDARD,
                 line_type=Line_Type.RL):
        if not isinstance(network, net.System):
            raise Exception("network must be an object of class Network of PyVolt")

        if not isinstance(measurement_set, measurement.MeasurementSet):
            raise Exception("measurement_set must be an object of class MeasurementSet of PyVolt")
        
        if not isinstance(critical_nodes, list):
            raise Exception("critical_nodes must be passed as a list object")
        
        if not all(isinstance(item, net.Node) for item in critical_nodes):
            raise Exception("all items of critical_nodes must be of PyVolt Network Node object type")
        
        if not isinstance(model_type, DpDse_Model):
            raise Exception("model_type must be an object of class DpDse_Model")

        if not isinstance(line_type, Line_Type):
            raise Exception("line_type must be an object of class Line_Type")

        self.model_type = model_type
        self.network = network
        self.measurement_set = measurement_set
        self.critical_nodes = critical_nodes
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
        self.u = np.array(
            [])  # list of control inputs in the form [vmag, vph, P, Q] or [vmag, vp, Iinj_mag, Iinj_ph] (TODO: second one needs new measurement type current injection in PyVolt)
        self.x_est = np.array([])  # list of values of estimated states
        self.x_pred = np.array([])  # list of values of predicted states
        self.P_est = np.array([])  # estimation covariance
        self.P_pred = np.array([])  # prediction covariance
        self.vg = np.array([])  # array of generator voltage phasors measurement object as tuple (v_mag, v_ph).
        self.sl = np.array([])  # array of critical load power injection measurement object as tuple (p, q)
        self.il = np.array([])  # array of critical load current injection measurement object as tuple (il_mag, il_ph)
        self.z = np.array([])  # array of measurement objects which will be used for correction step
        self.Ac = np.array([]) # critical load incidence matrix size l by c, c is number of critical loads
        self.Ac_d = np.array([]) # non-critical load incidence matrix
        self.num_c = 0 # number of critical loads
        self.num_nc = 0 # number of non-critical load
    

    def get_critical_non_critical_index(self):
        get_EC_node_index = [load_node.index for load_node in self.network.get_EC_nodes()]
        critical_nodes_index = [c_node.index for c_node in self.critical_nodes]
        non_critical_nodes_index = [item for item in get_EC_node_index if item not in critical_nodes_index]

        return critical_nodes_index, non_critical_nodes_index
    
    def get_non_critical_nodes(self):
        all_load_nodes = self.network.get_EC_nodes()
        return [item for item in all_load_nodes if item not in self.critical_nodes]
    
    def create_load_incidence_matrix(self):
        get_EC_node_index = [load_node.index for load_node in self.network.get_EC_nodes()]
        num_nodes = self.num_g + self.num_l
        A = np.zeros((num_nodes, num_nodes))
        c_nodes_idx, nc_nodes_idx = self.get_critical_non_critical_index()

        for i in range(num_nodes):
            for j in range(num_nodes):
                if  j in c_nodes_idx and i == j:
                    A[i, j] = 1
                if j in nc_nodes_idx and i == j:
                    A[i, j] = 1

        self.Ac = A[np.ix_(get_EC_node_index, c_nodes_idx)]
        self.Ac_d = A[np.ix_(get_EC_node_index, nc_nodes_idx)]



    def initialize_dse(self):
        fo = 50  # TODO: where can we get network nominal frequency information
        w_o = 2 * np.pi * fo

        # set type of line and if load resistance is to be considered # TODO: How and where to specify load resistance?
        load_resistance = False

        # TODO: How to set P_rLoad and R_L??
        # Set base quantities
        Vbase = 12.66  # line-line voltage #TODO: How to set vbase?
        base_apparent_power = 1  # MVA
        Ibase = base_apparent_power / (np.sqrt(3) * Vbase)

        get_ES_node_index = [gen_node.index for gen_node in self.network.get_ES_nodes()]
        get_EC_node_index = [load_node.index for load_node in self.network.get_EC_nodes()]

        num_gen = len(get_ES_node_index)
        num_load = len(get_EC_node_index)
        num_nodes = num_gen + num_load
        num_branch = self.network.get_branch_num()

        self.num_b = num_branch
        self.num_g = num_gen
        self.num_l = num_load
        self.num_c = len(self.critical_nodes)
        self.num_nc = num_load - self.num_c

        
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

        # Extracting R, L, C
        branch_r = [1e-15 if x == 0.0 else x for x in self.network.get_branch_R()]
        cables_R = np.diag(branch_r)
        branch_x = [1e-15 if x == 0.0 else x for x in
                    self.network.get_branch_X()]
        cables_L = np.diag(branch_x) / w_o
        cables_C = np.diag(np.dot(abs(PS_A_L), self.network.get_branch_BCH()) * 0.5) / w_o

        # create critical and non critical load incidence matrix Ac and Ac_d respectively
        self.create_load_incidence_matrix()

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
            A31 = np.dot(-np.eye(num_load), PS_A_L)
        A32 = np.zeros((num_load, num_branch))
        if self.line_type == Line_Type.PI and load_resistance is True:
            A33 = -np.linalg.inv(np.dot(cables_C, R_L))
        elif load_resistance is False:
            A33 = np.zeros((num_load, num_load))  # when load resistance is not to be considered
        elif self.line_type == Line_Type.RL and load_resistance is True:
            A33 = -np.linalg.inv(np.dot(np.eye(num_load), R_L))
        if self.line_type == Line_Type.PI:
            A34 = w_o * np.eye(num_load)
        else:
            A34 = np.zeros((num_load, num_load))  # when shunt capacitance is not to be considered

        A41 = np.zeros((num_load, num_branch))
        if self.line_type == Line_Type.PI:
            A42 = np.dot(-np.linalg.inv(cables_C), PS_A_L)
        else:
            A42 = np.dot(-np.eye(num_load), PS_A_L)
        if self.line_type == Line_Type.PI:
            A43 = -w_o * np.eye(num_load)
        else:
            A43 = np.zeros((num_load, num_load))  # when shunt capacitance is not to be considered
        if self.line_type == Line_Type.PI and load_resistance is True:
            A44 = -np.linalg.inv(np.dot(cables_C, R_L))
        elif load_resistance is False:
            A44 = np.zeros((num_load, num_load))  # when load resistance is not to be considered
        elif self.line_type == Line_Type.RL and load_resistance is True:
            A44 = -np.linalg.inv(np.dot(np.eye(num_load), R_L))

        # Creating B matrix
        B11 = np.dot(np.linalg.inv(cables_L), PS_A_G.transpose())
        B12 = np.zeros((num_branch, num_gen))
        B13 = np.zeros((num_branch, self.num_c))
        B14 = np.zeros((num_branch, self.num_c))

        B21 = np.zeros((num_branch, num_gen))
        B22 = np.dot(np.linalg.inv(cables_L), PS_A_G.transpose())
        B23 = np.zeros((num_branch, self.num_c))
        B24 = np.zeros((num_branch, self.num_c))

        B31 = np.zeros((num_load, num_gen))
        B32 = np.zeros((num_load, num_gen))
        if self.line_type == Line_Type.PI:
            B33 = -np.linalg.inv(cables_C) @ self.Ac
        else:
            B33 = -np.eye(num_load) @ self.Ac
        B34 = np.zeros((num_load, self.num_c))

        B41 = np.zeros((num_load, num_gen))
        B42 = np.zeros((num_load, num_gen))
        B43 = np.zeros((num_load, self.num_c))
        if self.line_type == Line_Type.PI:
            B44 = -np.linalg.inv(cables_C) @ self.Ac
        else:
            B44 = -np.eye(num_load) @ self.Ac

        # unknown control input matrix

        G11 = np.zeros((num_branch, self.num_nc))
        G12 = np.zeros((num_branch, self.num_nc))
        G21 = np.zeros((num_branch, self.num_nc))
        G22 = np.zeros((num_branch, self.num_nc))
        if self.line_type == Line_Type.PI:
            G31 = -np.linalg.inv(cables_C) @ self.Ac_d
        else:
            G31 = -np.eye(num_load) @ self.Ac_d
        G32 = np.zeros((num_load, self.num_nc))
        G41 = np.zeros((num_load, self.num_nc))
        if self.line_type == Line_Type.PI:
            G42 = -np.linalg.inv(cables_C) @ self.Ac_d
        else:
            G42 = -np.eye(num_load) @ self.Ac_d

        SS_A = np.bmat([[A11, A12, A13, A14], [A21, A22, A23, A24], [A31, A32, A33, A34], [A41, A42, A43, A44]])
        SS_B = np.bmat([[B11, B12, B13, B14], [B21, B22, B23, B24], [B31, B32, B33, B34], [B41, B42, B43, B44]])
        SS_G = np.bmat([[G11, G12], [G21, G22], [G31, G32], [G41, G42]])
        SS_Bcomb = np.bmat([[B11, B12, B13, B14, G11, G12], [B21, B22, B23, B24, G21, G22], [B31, B32, B33, B34, G31, G32], [B41, B42, B43, B44, G41, G42]])
        SS_A_disc = spy.linalg.expm(SS_A * self.time_step)
        SS_Bcomb_disc = np.dot(
            np.dot(np.linalg.inv(SS_A), (spy.linalg.expm(np.dot(SS_A, self.time_step)) - np.eye(np.shape(SS_A)[0]))),
            SS_Bcomb)
        SS_B_disc = SS_Bcomb_disc[:, :self.num_g + self.num_c]
        SS_G_disc = SS_Bcomb_disc[:, -self.num_nc:]


        if self.model_type == DpDse_Model.STANDARD:
            # create A and B matrix for standard case, i.e.,
            # SV: branch currents and load voltages in rectangular form
            # u : generator voltages and critical load current injections in rectangular form
            # uu : unknown control inputs - non-critical load current injections in rectangular form
            # either load current measurements are available or power injection measurements are available
            self.num_sv = 2 * num_branch + 2 * num_load
            self.num_u = 2 * num_gen + 2 * self.num_c 
            self.num_uu = 2*self.num_nc # number of unknown inputs 
            self.Act = SS_A
            self.Bct = SS_B
            self.Gct = SS_G
            self.Adt = SS_A_disc
            self.Bdt = SS_B_disc
            self.Gdt = SS_G_disc

        # TODO: implement for other DSE model augmented and load power
        # extract the control variables into self.vg and self.il or self.sl private variables
        self.separate_inputs()


        # initialize SV using static se or power flow
        try:
            # TODO: the current injection as measurement is not considered in PyVolt! Pyvolt will throw error!
            static_se_results = nv_state_estimator.DsseCall(self.network, self.measurement_set)
            self.initialize_sv(static_se_results)
        except:
            print("PyVolt SE threw error! State variables are initialized from power flow results!")
            static_se_results, num_iter = nv_powerflow.solve(self.network)
            self.initialize_sv(static_se_results)
        
        # initialize initial estimation covariance
        self.P_est = 1e-10 * np.ones((self.num_sv, self.num_sv))

    def check_ss_consistency(self):
        # TODO: This consistency check is only valid for STANDARD model type
        if self.model_type == DpDse_Model.STANDARD:
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

            c_idx, nc_idx = self.get_critical_non_critical_index()
            pf_il_c = np.array(pf_node_currents)[c_idx]
            pf_il_nc = np.array(pf_node_currents)[nc_idx]
            pf_ibr = np.array(pf_br_currents)
            # converting line-line voltage to phase voltage: the state-space is built for per phase calculation
            x_init_pf_conv = np.concatenate(
                (pf_ibr.real, pf_ibr.imag, pf_vl.real / np.sqrt(3), pf_vl.imag / np.sqrt(3)))
            u_init_pf_conv = np.concatenate((pf_vg.real / np.sqrt(3), pf_vg.imag / np.sqrt(3), pf_il_c.real,
                                                  pf_il_c.imag, pf_il_nc.real, pf_il_nc.imag))
            # negated load currents as current injections from loads are considered and rearranging critical and then non-critical nodes
            u_init_pf_curr_inj = np.concatenate((pf_vg.real / np.sqrt(3), pf_vg.imag / np.sqrt(3), - pf_il_c.real,
                                                 - pf_il_c.imag, - pf_il_nc.real, - pf_il_nc.imag))

            combine_B = np.bmat([[self.Bct, self.Gct]])
            bu = np.dot(combine_B, u_init_pf_conv)
            # solve for x = inv(A)*Bu
            bu_curr_inj = np.dot(combine_B, u_init_pf_curr_inj)
            x_init_from_ss = np.dot(np.linalg.inv(self.Act), -bu_curr_inj.T)

            print("x_init direct pf: ", x_init_pf_conv)
            print("x_init_state space: ", x_init_from_ss)

            print("Ax: ", np.dot(self.Act, x_init_pf_conv))
            print("Bu: ", bu)

            print("Ax_state space: ", np.dot(self.Act, x_init_from_ss).T)
            print("bu_curr_inj: ", bu_curr_inj)
        else:
            print("consistency check is implemented only for STANDARD type")

    def separate_inputs(self):
        # from the entire measurement_set, extract and separate the measurements which form control variables
        # (vg, il or sl) and rest of the measurements which form observations used for correct step

        v_mag = []  # list to extract generator voltage magnitude measurement objects
        v_ph = []  # list to extract generator voltage phase measurement objects
        p = []  # list to extract load real power injection measurement objects
        q = []  # list to extract load reactive power injection measurement objects
        i_mag = []  # list to extract load current injection magnitude measurement objects
        i_ph = []  # list to extract load current injection phase measurement objects
        z = []  # list of rest of the measurement objects for use in correction step
        z_val = []
        gen_uuid = [gen_node.uuid for gen_node in self.network.get_ES_nodes()]
        load_uuid = [load_node.uuid for load_node in self.network.get_EC_nodes()]

        # first extract all measurements of a node of type voltage and power injections of generator nodes and load
        # nodes respectively, the remaining goes into z list
        for meas in self.measurement_set.measurements:
            if meas.element_type == measurement.ElemType.Node:
                if meas.meas_type == measurement.MeasType.Vpmu_mag and meas.element.uuid in gen_uuid:
                    v_mag.append(meas)
                if meas.meas_type == measurement.MeasType.Vpmu_phase and meas.element.uuid in gen_uuid:
                    v_ph.append(meas)
                if meas.meas_type == measurement.MeasType.Sinj_real and meas.element.uuid in load_uuid:
                    p.append(meas)
                if meas.meas_type == measurement.MeasType.Sinj_imag and meas.element.uuid in load_uuid:
                    q.append(meas)
                if meas.meas_type == measurement.MeasType.Ipmu_inj_mag and meas.element.uuid in load_uuid:
                    i_mag.append(meas)
                if meas.meas_type == measurement.MeasType.Ipmu_inj_phase and meas.element.uuid in load_uuid:
                    i_ph.append(meas)
        
        all_u = v_mag + v_ph + p + q + i_mag + i_ph
        z = [item for item in self.measurement_set.measurements if item not in all_u]

        z_vals = [i.meas_value for i in z]
        print("vals in z_new: ", z_vals)

        # Create a dictionary to map uuid to measurement objects
        vmag_dict = {m.element.uuid: m for m in v_mag}
        vph_dict = {m.element.uuid: m for m in v_ph}
        p_dict = {m.element.uuid: m for m in p}
        q_dict = {m.element.uuid: m for m in q}
        i_inj_mag_dict = {m.element.uuid: m for m in i_mag}
        i_inj_ph_dict = {m.element.uuid: m for m in i_ph}

        # Rearrange gen voltages and according to the order of generator order (using uuid), same for load powers
        # list of control variables - (voltage magnitude, voltage phase angles),
        #                             (load real power and load imaginary power injections)
        if vmag_dict and vph_dict:
            self.vg = [(vmag_dict[x], vph_dict[x]) for x in gen_uuid]
        if p_dict and q_dict: # here it is assumed that real and reactive power are obtained together
            self.sl = [(p_dict[x.uuid], q_dict[x.uuid]) for x in self.critical_nodes] # only critical loads injection for known control inputs
        if i_inj_mag_dict and i_inj_ph_dict:
            self.il = [(i_inj_mag_dict[x.uuid], i_inj_ph_dict[x.uuid]) for x in self.critical_nodes] # only critical loads injection for known control inputs
        self.z = np.array(z)

        # TODO: check for power injections and current injections, the sign needs to be changed or not

    def assemble_inputs(self):
        # this function constructs u vector depending upon the kind of model in the required form
        # if STANDARD: [vg_re, vg_im, il_re, il_im]
        # if INPUT_LOAD_POWER: [vg_re, vg_im, il_calc_re, il_calc_im]
        # if AUGMENTED: [0.....0]
        #TODO: yet to check if this is working correctly!!
        vg_mag = np.array([m[0].meas_value for m in self.vg])
        vg_ph = np.array([m[1].meas_value for m in self.vg])
        vg_cmplx = phasor_complex(vg_mag, vg_ph)
        vg = np.concatenate((vg_cmplx.real, vg_cmplx.imag))
        if self.model_type == DpDse_Model.STANDARD:
            il_mag = np.array([m[0].meas_value for m in self.il])
            il_ph = np.array([m[1].meas_value for m in self.il])
            il_cmplx = phasor_complex(il_mag, il_ph)
            il = np.concatenate((il_cmplx.real, il_cmplx.imag))
            self.u = np.concatenate((vg, il))
        elif self.model_type == DpDse_Model.INPUT_LOAD_POWER:
            il_calc = self.power_to_current_injection()
            self.u = np.concatenate((vg, il_calc))
        elif self.model_type == DpDse_Model.AUGMENTED:
            self.u = np.zeros((self.num_sv, 1))  # no control variables, it is only X_dot = Ax

    def power_to_current_injection(self):
        # compute current injections from power injection measurements and estimated voltages from previous step
        #TODO: need to check if this function works!
        v_est = self.extract_load_voltages_estimation()
        v_est_re = v_est[:self.num_l]
        v_est_im = v_est[self.num_l:]
        p_inj = np.array([m[0].meas_value for m in self.sl])
        q_inj = np.array([m[1].meas_value for m in self.sl])
        v_sq = v_est_re * v_est_re + v_est_im * v_est_im
        il_re = (v_est_re * p_inj + v_est_im * q_inj) / v_sq
        il_im = (v_est_im * p_inj - v_est_re * q_inj) / v_sq
        il = np.concatenate((il_re, il_im))
        return il

    def initialize_sv(self, static_se_results):
        # Perform state estimation
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
        self.x_est = np.concatenate((ib_re, ib_im, vl_re, vl_im))

    def predict(self):
        # TODO: yet to check if this works correctly!!
        # assemble u vector depending upon the model type into required form
        self.assemble_inputs()
        print("self.u: ", self.u)
        for i in self.z:
            print(i.meas_value)
        # predict the states for next time-step
        # self.x_pred = self.Adt * self.x_est + self.Bdt * self.u.T
        self.x_pred = self.Adt @ self.x_est
        self.P_pred = self.Adt @ self.P_est @ (self.Adt).T

    
    def correct(self):
        # Step 1: look if new measurements received
        # Step 2: update the measurements vector
        # Step 3: build the measurement functions
        # Step 4: build the Jacobian
        # Step 5: build F matrix
        # Step 6: build R matrix
        # Step 7: build R_d matrix
        # Step 8: build M matrix
        # Step 9: build d vector of unknown inputs
        # Step 10: update prediction x_pred
        # Step 11: build Kalman gain
        # Step 12: estimate x_est

        
        return 1
        

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

    def hx_Ipmu_b_re(self):
        # measurement function for measurements of branch currents - real
        return self.extract_branch_currents_prediction()[:self.num_b]

    def hx_Ipmu_b_im(self):
        # measurement function for measurements of branch currents - imag
        return self.extract_branch_currents_prediction()[self.num_b: 2 * self.num_b]

    def H_Ipmu_b_re(self):
        # Jacobian for pmu measurements of branch currents
        print('To be implemented')

    def H_Ipmu_b_im(self):
        # Jacobian for pmu measurements of branch currents
        print('To be implemented')

    def hx_Vpmu(self):
        # meausrement function for pmu measurements of load voltages
        print('To be implemented')

    def H_Vpmu(self):
        # Jacobian for pmu measurements of load voltages
        print('To be implemented')
