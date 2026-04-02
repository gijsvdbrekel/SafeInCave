# Copyright 2025 The safeincave community.
#
# This file is part of safeincave.
#
# Licensed under the GNU GENERAL PUBLIC LICENSE, Version 3 (the "License"); you may not
# use this file except in compliance with the License.  You may obtain a copy
# of the License at
#
#     https://spdx.org/licenses/GPL-3.0-or-later.html
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
# WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  See the
# License for the specific language governing permissions and limitations under
# the License.
from abc import ABC, abstractmethod
import CoolProp.CoolProp as CP
import numpy as np


class Cavern(ABC):
    def __init__(self,
                 cavern_name: str,
                 fluid: str,
                 h_conv: float = None):
        self.cavern_name = cavern_name
        self.fluid = fluid
        self.h_conv = h_conv
        self.type = None

    def check_fluid(self) -> None:
        try:
            CP.AbstractState("HEOS", self.fluid)
        except ValueError:
            raise ValueError(f"Fluid '{self.fluid}' not recognized by CoolProp.")
        
    @abstractmethod
    def update_cavern(self, t: float) -> None:
        pass



class CavernHandler:
    def __init__(self):
        self.caverns_T = []
        self.caverns_PT = []

    def add_cavern(self, cavern: Cavern) -> None:
        if cavern.type == "Cavern_T":
            self.caverns_T.append(cavern)
        elif cavern.type == "Cavern_PT":
            self.caverns_PT.append(cavern)
        else:
            raise ValueError(f"Cavern type {cavern.type} not supported")

    def update_caverns(self, t: float):
        for cavern in self.caverns_T:
            cavern.update_cavern(t)

        for cavern in self.caverns_PT:
            cavern.update_cavern(t)


class Cavern_T(Cavern):
    def __init__(self, cavern_name: str, T_values: list, time_values: list, h_conv: float = None):
        super().__init__(cavern_name, None, h_conv)
        self.type = "Cavern_T"
        self.T_values = T_values
        self.time_values = time_values
        self.T = self.T_values[0]

    def update_cavern(self, t: float) -> None:
        self.T = np.interp(t, self.time_values, self.T_values)
        if self.T <= 0.0:
            raise ValueError(f"T must be > 0, got {self.T}")

class Cavern_PT(Cavern):
    def __init__(self, 
                 cavern_name: str,
                 fluid: str,
                 P_values: list,        # Gauge pressure values (Pa)
                 T_values: list,        # Temperature values (K)
                 time_values: list,
                 ref_pos: float = 0.0,
                 direction: int = 0,
                 g: float = -9.81,
                 h_conv: float = None,
                 P_atm: float = 101325.0 # Atmospheric pressure in Pa
                ):
        super().__init__(cavern_name, fluid, h_conv)
        self.type = "Cavern_PT"
        self.P_values = P_values
        self.T_values = T_values
        self.time_values = time_values
        self.P_atm = P_atm
        self.ref_pos = ref_pos
        self.direction = direction
        self.gravity = g
        self.__initialize()

    def __initialize(self) -> None:
        self.AS = CP.AbstractState("HEOS", self.fluid)
        self.P = self.P_values[0]
        self.T = self.T_values[0]
        self.AS.update(CP.PT_INPUTS, self.P, self.T)
        self.density = self.AS.rhomass()

    def update_cavern(self, t: float) -> None:
        self.P = np.interp(t, self.time_values, self.P_values)
        self.T = np.interp(t, self.time_values, self.T_values)
        if self.P <= 0.0:
            raise ValueError(f"P must be > 0, got {self.P}")
        if self.T <= 0.0:
            raise ValueError(f"T must be > 0, got {self.T}")
        self.AS.update(CP.PT_INPUTS, self.P + self.P_atm, self.T)
        self.density = self.AS.rhomass()



class Cavern_MassFlux(Cavern):
    def __init__(self, 
                 cavern_name: str,
                 fluid: str,
                 P_init: float,        # Initial gauge pressure (Pa)
                 T_init: float,
                 Mflux_values: list,
                 time_values: list,
                 ref_pos: float = 0.0,
                 direction: int = 0,
                 g: float = -9.81,
                 h_conv: float = None,
                 P_atm: float = 101325.0 # Atmospheric pressure in Pa
                 ):
        super().__init__(cavern_name, fluid, h_conv)
        self.type = "Cavern_MassFlux"
        self.Mflux_values = Mflux_values
        self.time_values = time_values
        self.Mflux = self.Mflux_values[0]
        self.ref_pos = ref_pos
        self.direction = direction
        self.gravity = g
        self.P_atm = P_atm
        self.P_init = P_init
        self.T_init = T_init
        self.__initialize()

    def __initialize(self) -> None:
        self.AS = CP.AbstractState("HEOS", self.fluid)
        self.AS.update(CP.PT_INPUTS, self.P_init + self.P_atm, self.T_init)
        self.density = self.AS.rhomass()

    def update_cavern(self, t: float) -> None:
        self.Mflux = np.interp(t, self.time_values, self.Mflux_values)
        if self.Mflux < 0.0:
            raise ValueError(f"Mass flux must be >= 0, got {self.Mflux}")

    