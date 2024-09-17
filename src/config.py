from enum import Enum

class DpDse_Mode(Enum):
    ONLINE = 1
    OFFLINE = 2

class DpDse_Output(Enum):
    STORE = 1
    STREAM = 2

class DpDse_Smoother(Enum):
    Yes = 1
    No = 2

class DSEConfig:
    def __init__(self, mode, output, duration, smoother):
        self.mode = mode  # DSE in offline or online mode
        self.output = output  # Output to be streamed or stored
        self.duration = duration  # duration for DSE to run; applicable in both modes
        self.smoother = smoother  # smoother enabled or not


    # Getter for mode
    @property
    def mode(self):
        return self._mode

    @property
    def output(self):
        return self._output

    @property
    def duration(self):
        return self._duration

    @property
    def smoother(self):
        return self._smoother


    # Setter for mode with validation
    @mode.setter
    def mode(self, value):
        if value not in [DpDse_Mode.ONLINE, DpDse_Mode.OFFLINE]:
            raise ValueError("Mode must be of object DpDse_Mode.")
        self._mode = value

    @output.setter
    def output(self, value):
        if value not in [DpDse_Output.STORE, DpDse_Output.STREAM]:
            raise ValueError("Output must be of object DpDse_Output.")
        self._output = value

    @duration.setter
    def duration(self, value):
        if value == 0.00:
            raise ValueError("Duration of DSE execution is set to 0. DSE will not be executed")
        self._duration = value

    @smoother.setter
    def smoother(self, value):
        if value not in [DpDse_Smoother.Yes, DpDse_Smoother.No]:
            raise ValueError("smoother must be of object DpDse_Smoother.")
        self._smoother = value



