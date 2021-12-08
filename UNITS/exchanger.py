# heat exchanger unit class handler
from stream import Stream
import copy
class Exchanger():
    """
    perform the heat balance and energy calculations on given streams
    """
    def __init__(self, name, hotFluid: Stream, coldFluid: Stream, approach_temp=5, flow="counter-current"):
        """
        name: name of the heat exchanger for identification purpose
        hotFluid: hot fluid of Stream Class object
        coldFluid: cold fluid of Stream Class object
        optional approach_temp: minimum temperature difference allowed for heat transfer
        optional flow: flow arrangement inside exchanger, "co-current" or "counter-current"

        RETURNS
        ----------------------------
        (hotDischarge, coldDischarge, exchanger_duty) tuple with updated tempertures
        * only implmenting it for mixture as this is our use case
        """
        self.name = name
        self.hotFluid = hotFluid
        self.coldFluid = coldFluid
        self.operation_mode = flow
        self.approach_temp = approach_temp
        (hot_discharge, cold_discharge, dh) = self.solve()
        self.hot_discharge = hot_discharge
        self.cold_discharge = cold_discharge
        self.duty = dh


    def solve(self):
        print("="*50)
        print("Exchanger feed:\n")
        print("Hot Fluid: \n")
        self.hotFluid.show()
        print("Cold Fluid: \n")
        self.coldFluid.show()


        # checking input for validation
        if (self.hotFluid.temp - self.coldFluid.temp) <= self.approach_temp:
            # heat transfer is not possible, ask user to change the inputs and rerun
            print("given conditions cannot satisfy minimum approach temperture")
            print("please update input streams or minimum approach temperture")
            raise ValueError("improper temperture specifications")
        else:
            # tempertures are valid

            if self.operation_mode == "counter-current":
                # assume hot side discharge temp
                hot_dischargeT = self.coldFluid.temp + self.approach_temp
                # assume cold side discharge temp
                cold_dischargeT = self.hotFluid.temp - self.approach_temp
            elif self.operation_mode == "co-current":
                # has to be done iteratively
                pass
            # calculate enthalpy change for hot and cold
            dh_hot = self.hotFluid.delta_h(hot_dischargeT)    # w
            dh_cold = self.coldFluid.delta_h(cold_dischargeT)  # w
            # compare the two changes and adjust temp so the enthalpy changes are equal to the smaller value
            if abs(dh_hot) < abs(dh_cold):
                # cold stream will not reach minimum approach temp
                cold_dischargeT = self.coldFluid.discharge_temp(-dh_hot)
            else:
                # hot stream will not meet approach conditions
                hot_dischargeT = self.hotFluid.discharge_temp(-dh_cold)

        hot_discharge = copy.deepcopy(self.hotFluid)
        hot_discharge.temp = hot_dischargeT
        cold_discharge = copy.deepcopy(self.coldFluid)
        cold_discharge.temp = cold_dischargeT
        print("\n----------------------------------------------------")
        print("Exchanger Discharge:\n")
        print("Hot Discharge Temperture: ", hot_dischargeT)
        print("Cold Discharge Temperture: ", cold_dischargeT)
        print("Exchanger Duty: ", abs(dh_hot))
        return (hot_discharge, cold_discharge, abs(dh_hot))

