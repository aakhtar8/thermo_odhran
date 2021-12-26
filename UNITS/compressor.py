# Compressor class to simulate compressors, discharge temperture and work required will be calculated
from stream import Stream

class Compressor:
    def __init__(self, name: str, stream: Stream, discharge_pres = None, pres_ratio = None, isentropic = True):
        """
        Compressor class will take a stream and perform isentropic or polytropic calculations for given compression specfications
        -------------------------------------
        PARAMETERS
        name: (str) name of the unit for identification
        stream: Stream class instance, specifying the fluid which compressor is need to process
        Optional discharge_pres(Pa): discharge pressure that is required of the compressor
        Optional pressure_ratio: ratio of discharge to suction pressure
        * Either pressure ratio of discharge pressure has to be specified
        OPtional isentropic (boolean): True if compression is to be assumed isentropic, false for polytropic compression
        ----------------------------------------
        ATTRIBUTES
        name: name of unit
        feed: stream class instance representing feed for the compressor
        discharge: stream class representing discharge condtions and fluid for compressor
        isentropic: True or False depending upon the simulation environment desired
        work: work required by the compressor for given discharge specification
        pres_ratio: ratio of discharge to suction pressure
        """
        self.name = name
        self.feed = stream
        if discharge_pres is None:  # check pressure ratio
            if pres_ratio is None:  # throw an error
                raise ValueError(f"Invalid discharge specification for compressor {self.name}")
            else:
                self.pres_ratio = pres_ratio
        else:
            if pres_ratio is not None:  # bothe are specified, system becomes overspecified
                raise ValueErrorf(f"Invalid discharge specification for compressor {self.name}, system is over-specified")
            else:
                self.pres_ratio = discharge_pres/self.feed.pres
        self.isentropic = isentropic
        self.discharge = pass
        self.work = pass
    

    def polytropic_temp(self):
        discharge_pres = self.feed.pres*self.pres_ratio
        density = 