# Throttling Valve class handler
from stream import Stream
import copy
class Throttling_Valve():
    """
    perform throttling simulation via PH flash in EOS modelling of the system
    """
    def __init__(self, stream: Stream, discharge_pres=None, pressure_reduction=None):
        """
        for given discharge pressure throttle valve will provide discharge temperture of the stream
        PARAMETERS
        stream: A Stream class object
        discharge_pres[OPTIONAL]: throttling valve discharge pressure in Pa
        pressure_reduction[OPTIONAL]: ratio of discharge to feed to pressure
        """
        self.feed = stream
        print("="*50)
        print("\nValve Feed: \n")
        stream.show()
        if pressure_reduction is not None:
            self.discharge_pres = stream.pres*pressure_reduction
        else:
            self.discharge_pres = discharge_pres
        pt_flash = stream.material.flash(P=stream.pres, T=stream.temp, zs = stream.composition)
        enthalpy = pt_flash.H()
        throttle_T = stream.material.flash(P=discharge_pres, H=enthalpy, zs = stream.composition).T
        discharge = copy.deepcopy(self.feed)
        discharge.temp =throttle_T
        discharge.pres = discharge_pres
        self.discharge = discharge
        print("\n----------------------------------------------------")
        print("Valve Discharge: \n")
        discharge.show()