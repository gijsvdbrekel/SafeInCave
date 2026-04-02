from safeincave.CavernBC import CavernHandler, Cavern_PT, Cavern_T
import CoolProp.CoolProp as CP
import torch as to
import numpy as np
import unittest

class Test_T(unittest.TestCase):
    def setUp(self):
        self.cavern_name = "Cave-1"
        self.cavern_1 = Cavern_T(self.cavern_name, [303, 400.0], [0.0, 10.0])

    def test_initialize(self):
        self.assertEqual(self.cavern_1.T, 303)


class Test_PT(unittest.TestCase):
    def setUp(self):
        self.fluid_1 = "Water"
        self.cave_name_1 = "Cavern_1"
        self.cave_1 = Cavern_PT(self.cave_name_1, self.fluid_1, [1e5, 1e5], [303, 400.0], [0.0, 10.0])
        
        self.fluid_2 = "Methane"
        self.cave_name_2 = "Cavern_2"
        self.cave_2 = Cavern_PT(self.cave_name_2, self.fluid_2, [1e5, 1e5], [303, 400.0], [0.0, 10.0])

        self.cavern_set = CavernHandler()
        self.cavern_set.add_cavern(self.cave_1)
        self.cavern_set.add_cavern(self.cave_2)

    def test_initialize(self):
        self.assertEqual(self.cave_1.AS.phase(), 0)  # CP.iphase_liquid
        self.assertEqual(self.cave_1.cavern_name, self.cave_name_1)
        self.assertEqual(self.cave_1.fluid, self.fluid_1)
        self.assertEqual(self.cave_1.AS.phase(), CP.get_phase_index("phase_liquid"))
        self.assertEqual(self.cave_1.density, 995.6940736483174)
        
        self.assertEqual(self.cave_2.AS.phase(), 2)  # CP.iphase_liquid
        self.assertEqual(self.cave_2.cavern_name, self.cave_name_2)
        self.assertEqual(self.cave_2.fluid, self.fluid_2)
        self.assertEqual(self.cave_2.AS.phase(), CP.get_phase_index("phase_supercritical_gas"))
        self.assertEqual(self.cave_2.density, 0.6378365378697821)

    def test_update_cavern(self):
        self.cave_1.update_cavern(5.0)
        self.assertEqual(self.cave_1.AS.phase(), CP.get_phase_index("phase_liquid"))
        self.assertEqual(self.cave_1.density, 972.8112951044953)
        
        self.cave_1.update_cavern(10.0)
        self.assertEqual(self.cave_1.AS.phase(), CP.get_phase_index("phase_gas"))
        self.assertEqual(self.cave_1.density, 0.5476054152259427)

    def test_cavern_handler(self):
        self.assertEqual(len(self.cavern_set.caverns_PT), 2)
        self.cavern_set.update_caverns(2.0)
        
        self.assertEqual(self.cavern_set.caverns_PT[0].AS.phase(), CP.get_phase_index("phase_liquid"))
        self.assertEqual(self.cavern_set.caverns_PT[0].density, 988.3718407745588)
        
        self.assertEqual(self.cavern_set.caverns_PT[1].AS.phase(), CP.get_phase_index("phase_supercritical_gas"))
        self.assertEqual(self.cavern_set.caverns_PT[1].density, 0.5992487366380541)
