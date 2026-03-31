
from abc import ABC
import CoolProp.CoolProp as CP
import numpy as np


class Cavern(ABC):
    def __init__(self):
        self.name = None

class Cavern_PT(Cavern):
    def __init__(self, name: str, fluid: str, P_values: list, T_values: list, time_values: list):
        self.name = name
        self.fluid = fluid
        self.P_values = P_values
        self.T_values = T_values
        self.time_values = time_values
        self.AS = CP.AbstractState("PR", fluid)
        self.initialize()

    def initialize(self):
        self.P = self.P_values[0]
        self.T = self.T_values[0]
        self.AS.update(CP.PT_INPUTS, self.P, self.T)
        self.rho = self.AS.rhomass()

    def update_cavern(self, t: float):
        self.P = np.interp(t, self.time_values, self.P_values)
        self.T = np.interp(t, self.time_values, self.T_values)
        self.AS.update(CP.PT_INPUTS, self.P, self.T)
        self.rho = self.AS.rhomass()
