from pyvolt import measurement
import pandas as pd
import numpy as np


# TODO: How to know if node current and node power has positive value for consumption (out-of grid) or injection (into-grid)???
# In PyVolt program, it is assumed that node current and node power has positive if it is consumption (out-of grid), therefore negated to make it injection - in assemble_u function of dpdse.py
# From RTDS, we get directly the injection currents, therefore not negated here. 
#TODO: Need a standard way to fix this.
def update_measurement(element_uuid, meas_type, meas_data, mapping, value_in_pu=True):
        """
        to update meas_value of a specific measurement object in the measurements array
        """
        load_power_input_type = [measurement.MeasType.Sinj_real, measurement.MeasType.Sinj_imag]
        load_current_input_type = [measurement.MeasType.Ipmu_inj_mag, measurement.MeasType.Ipmu_inj_phase]

        key = (element_uuid, meas_type)
        if key in mapping:
            meas = mapping[key]
        # only update measurements that are already included in the measurements set
            if meas.element.uuid == element_uuid:
                # special case for voltage magnitude: SOGNO interface only knows Vpmu_mag while measurement set distincts between Vpmu_mag and V_mag
                if meas_type == measurement.MeasType.Vpmu_mag and (meas.meas_type == measurement.MeasType.Vpmu_mag or meas.meas_type == measurement.MeasType.V_mag):
                    # volt pu conversion assuming that meas_data from device are in volts and single-phase value according to sogno interface
                    # while baseVoltage from CIM in [kV] and three-phase value
                    if not value_in_pu:
                        meas_value_pu = meas_data / (meas.element.baseVoltage / np.sqrt(3) )
                        meas_value_act = meas_data
                    else:
                        meas_value_pu = meas_data
                        meas_value_act = meas_data * (meas.element.baseVoltage / np.sqrt(3) )
                    print("Updating measurement value for {} of type {} from {} to {}, value_in_pu {} ".format(meas.element.uuid, str(meas.meas_type), meas.meas_value_act, meas_value_act, value_in_pu))                    
                    meas.meas_value = meas_value_pu 
                    meas.meas_value_act = meas_value_act                
                # special case for voltage magnitude: SOGNO interface only knows Ipmu_mag while measurement set distincts between Ipmu_mag and I_mag
                elif (meas_type == measurement.MeasType.Ipmu_mag or meas_type == measurement.MeasType.Ipmu_inj_mag) and (meas.meas_type == measurement.MeasType.Ipmu_mag or meas.meas_type == measurement.MeasType.I_mag or meas.meas_type == measurement.MeasType.Ipmu_inj_mag):
                    # current pu conversion assuming that meas_data from device are in A and single-phase value according to sogno interface
                    # while base_current in [kA]
                    if not value_in_pu:
                        meas_value_pu = meas_data / (meas.element.base_current )
                        meas_value_act = meas_data
                    else:
                        meas_value_pu = meas_data
                        meas_value_act = meas_data * (meas.element.base_current )
                    print("Updating measurement value for {} of type {} from {} to {}, value_in_pu {}".format(meas.element.uuid, str(meas.meas_type), meas.meas_value_act, meas_value_act, value_in_pu))
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
                    print("Updating measurement value for {} of type {} from {} to {}, value_in_pu {}".format(meas.element.uuid, str(meas.meas_type), meas.meas_value_act, meas_value_act, value_in_pu))
                    meas.meas_value = meas_value_pu
                    meas.meas_value_act = meas_value_act
                elif (meas_type == meas.meas_type and (meas_type == measurement.MeasType.Vpmu_phase or meas_type == measurement.MeasType.Ipmu_phase or meas_type == measurement.MeasType.Ipmu_inj_phase)):
                    print("Updating measurement value for {} of type {} from {} to {}, value_in_pu {}".format(meas.element.uuid, str(meas.meas_type), meas.meas_value_act, meas_data, value_in_pu))
                    if meas_type == measurement.MeasType.Ipmu_inj_phase:  
                        meas.meas_value = meas_data    # for negating for current injection
                        meas.meas_value_act = meas_data 
                    else:
                        meas.meas_value = meas_data
                        meas.meas_value_act = meas_data
        else:
            print(f"with uuid={element_uuid} and type={meas_type} not found")

