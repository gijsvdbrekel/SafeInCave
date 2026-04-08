from safeincave.CavernBC import CavernHandler, Cavern_PT, Cavern_T, CavernVolumeComputer
from safeincave.Utils import create_field_nodes
import CoolProp.CoolProp as CP
from safeincave import GridHandlerGMSH
import dolfinx as do
import numpy as np
import os
import unittest


class Test_CavernVolumeComputer(unittest.TestCase):
    def setUp(self):
        self.grid = GridHandlerGMSH("geom", os.path.join("files", "cube_caverns"))

        self.__calculate_expected_volume()
        self.__expected_normals()

    def __calculate_expected_volume(self):
        Lz = 1.0
        Ly = 2.0
        Lx = 1.0
        height = Lz/2
        R = 0.1*(Lx+Ly)/2
        self.volume_expected = 2*R*2*R*height

    def __expected_normals(self):
        self.normals_full = np.array([
                [ 0.00000000e+00, -1.00000000e+00,  0.00000000e+00],
                [ 0.00000000e+00, -1.00000000e+00,  0.00000000e+00],
                [ 0.00000000e+00,  0.00000000e+00,  1.00000000e+00],
                [ 0.00000000e+00,  0.00000000e+00,  1.00000000e+00],
                [-1.00000000e+00, -7.40148699e-16,  0.00000000e+00],
                [-1.00000000e+00,  0.00000000e+00,  4.44089212e-16],
                [-0.00000000e+00, -1.00000000e+00,  0.00000000e+00],
                [ 0.00000000e+00, -1.00000000e+00,  0.00000000e+00],
                [ 0.00000000e+00,  0.00000000e+00,  1.00000000e+00],
                [ 1.00000000e+00,  0.00000000e+00, -0.00000000e+00],
                [ 1.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                [ 0.00000000e+00,  0.00000000e+00,  1.00000000e+00],
                [-0.00000000e+00,  0.00000000e+00, -1.00000000e+00],
                [-1.00000000e+00, -0.00000000e+00, -4.44089208e-16],
                [ 0.00000000e+00, -0.00000000e+00, -1.00000000e+00],
                [ 1.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                [ 0.00000000e+00,  0.00000000e+00, -1.00000000e+00],
                [-1.00000000e+00,  7.40148667e-16,  0.00000000e+00],
                [ 0.00000000e+00,  1.00000000e+00, -4.44089213e-16],
                [ 7.40148683e-16,  1.00000000e+00, -0.00000000e+00],
                [ 1.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                [-7.40148683e-16,  1.00000000e+00,  0.00000000e+00],
                [ 0.00000000e+00, -0.00000000e+00, -1.00000000e+00],
                [ 0.00000000e+00,  1.00000000e+00,  4.44089207e-16],
        ])
        self.normals_half = np.array([
                [-0.,  0.,  1.],
                [ 0.,  0.,  1.],
                [ 1.,  0., -0.],
                [ 1.,  0.,  0.],
                [-0.,  0.,  1.],
                [ 0.,  0.,  1.],
                [-1.,  0.,  0.],
                [-1.,  0.,  0.],
                [ 0.,  0., -1.],
                [-0.,  0., -1.],
                [ 0.,  0., -1.],
                [ 0., -0., -1.],
                [ 0.,  1.,  0.],
                [ 0.,  1., -0.],
                [ 0.,  1.,  0.],
                [ 0.,  1.,  0.]
        ])
        self.normals_quarter = np.array([
                [ 0.,  0.,  1.],
                [ 0.,  0.,  1.],
                [ 0., -0.,  1.],
                [-0.,  0.,  1.],
                [ 0., -1.,  0.],
                [ 0., -1.,  0.],
                [-1.,  0.,  0.],
                [-1.,  0.,  0.],
                [-0.,  0., -1.],
                [ 0., -0., -1.],
                [ 0.,  0., -1.],
                [ 0.,  0., -1.]
        ])

    def test_compute_volume(self):
        cvc = CavernVolumeComputer(self.grid, "Cavern_full")
        volume_full = cvc.compute()
        self.assertAlmostEqual(volume_full, self.volume_expected, delta=1e-8)
        
        cvc_half = CavernVolumeComputer(self.grid, "Cavern_half", [10.0, 2.0, -780.0], sym_scale=2)
        volume_half = cvc_half.compute()
        self.assertAlmostEqual(volume_half, self.volume_expected, delta=1e-8)
        
        cvc_quarter = CavernVolumeComputer(self.grid, "Cavern_quarter", [0.0, 0.0, 0.5], sym_scale=4)
        volume_quarter = cvc_quarter.compute()
        self.assertAlmostEqual(volume_quarter, self.volume_expected, delta=1e-8)

    def test_normals(self):
        cvc = CavernVolumeComputer(self.grid, "Cavern_full")
        normals = cvc.calculate_normals()
        np.testing.assert_allclose(normals, self.normals_full, rtol=1e-3)
        
        cvc = CavernVolumeComputer(self.grid, "Cavern_half")
        normals = cvc.calculate_normals()
        np.testing.assert_allclose(normals, self.normals_half, rtol=1e-3)
        
        cvc = CavernVolumeComputer(self.grid, "Cavern_quarter")
        normals = cvc.calculate_normals()
        np.testing.assert_allclose(normals, self.normals_quarter, rtol=1e-3)

    def test_compute_volume_with_displacement(self):
        cvc = CavernVolumeComputer(self.grid, "Cavern_full")
        CG1_3x1 = do.fem.functionspace(cvc.grid.mesh, ("Lagrange", 1, (self.grid.domain_dim, )))
        u = do.fem.Function(CG1_3x1)
        u_array = u.x.array.reshape(-1, 3).copy()

        fun_ux = lambda x, y, z: 0.5 - x*0.2
        u_array[:, 0] = create_field_nodes(cvc.grid, fun_ux)
        u.x.array[:] = u_array.flatten()
        volume_computed = cvc.compute(u)

        R = 0.15
        p1 = 0.5 - R
        p2 = 0.5 + R
        u1 = 0.5 - p1*0.2
        u2 = 0.5 - p2*0.2
        D = (p2 + u2) - (p1 + u1)
        height = 0.75 - 0.25
        V_expected = D*2*R*height
        
        self.assertAlmostEqual(volume_computed, V_expected, delta=1e-8)



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
