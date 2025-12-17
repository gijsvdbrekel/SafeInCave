# import os
# import sys
# sys.path.append(os.path.join("..", "..", "..", "safeincave"))
# from Grid2 import GridHandlerGMSH
# from mpi4py import MPI
# import ufl
# import dolfinx as do
# import torch as to
# import numpy as np
# from petsc4py import PETSc
# import Utils as utils
# from MaterialProps import *
# from MomentumEquation import LinearMomentum
# import MomentumBC as momBC
# from OutputHandler import SaveFields
# from Simulators import Simulator_M
# from TimeHandler import TimeController, TimeControllerParabolic
# import time

from safeincave.Utils import GPa, MPa, day, hour, create_field_elems, create_field_nodes
import safeincave as sf
import safeincave.Utils as ut
import safeincave.MomentumBC as momBC
from mpi4py import MPI
import dolfinx as do
import os
import sys
import ufl
import torch as to
import numpy as np
from petsc4py import PETSc
import time

GPa = ut.GPa
MPa = ut.MPa
day = ut.day


def get_geometry_parameters(path_to_grid):
	f = open(os.path.join(path_to_grid, "geom.geo"), "r")
	data = f.readlines()
	ovb_thickness = float(data[8][len("ovb_width = "):-2])
	salt_thickness = float(data[9][len("salt_width = "):-2])
	hanging_wall = float(data[10][len("roof_width = "):-2])
	return ovb_thickness, salt_thickness, hanging_wall


def create_input_file(P_min, P_max, rho_ob, E_ob, nu_ob, rho_salt, E_salt, nu_salt, A_salt, n_salt, case_folder):

	def main():
		comm = MPI.COMM_WORLD
		comm.Barrier()
		if MPI.COMM_WORLD.rank == 0:
		    start_time = MPI.Wtime()

		# Read grid
		grid_path = os.path.join("grid")
		grid = sf.GridHandlerGMSH("geom", grid_path)

		# Define output folder
		output_folder = os.path.join("output", case_folder)

		# Define momentum equation
		mom_eq = sf.LinearMomentum(grid, theta=0.0)

		# Define solver
		mom_solver = PETSc.KSP().create(grid.mesh.comm)
		mom_solver.setType("bicg")
		mom_solver.getPC().setType("asm")
		mom_solver.setTolerances(rtol=1e-12, max_it=100)
		mom_eq.set_solver(mom_solver)

		# Define material properties
		mat = sf.Material(mom_eq.n_elems)

		# Extract region indices
		ind_salt = grid.region_indices["SALT"]
		ind_ovb = grid.region_indices["OVERBURDEN"]

		# Set material density
		salt_density = rho_salt
		ovb_density = rho_ob
		gas_density = 0.083
		R = 8.32
		rho = to.zeros(mom_eq.n_elems, dtype=to.float64)
		rho[ind_salt] = salt_density
		rho[ind_ovb] = ovb_density
		mat.set_density(rho)

		# Constitutive model
		E0 = to.zeros(mom_eq.n_elems)
		E0[ind_salt] = E_salt #102*GPa
		E0[ind_ovb] = E_ob #180*GPa
		nu0 = to.zeros(mom_eq.n_elems)
		nu0[ind_salt] = nu_salt
		nu0[ind_ovb] = nu_ob
		spring_0 = sf.Spring(E0, nu0, "spring")

		# Create Kelvin-Voigt viscoelastic element
		eta = 105e11*to.ones(mom_eq.n_elems)
		E1 = 10*GPa*to.ones(mom_eq.n_elems) 
		nu1 = 0.3*to.ones(mom_eq.n_elems)
		kelvin = sf.Viscoelastic(eta, E1, nu1, "kelvin")


		# Constant Temperature
		def T_field_fun(x, y, z):
		    return 300.0
		T0_field = create_field_elems(grid, T_field_fun)
		mom_eq.set_T0(T0_field)
		mom_eq.set_T(T0_field)

		# Create creep
		Q = 51600*to.ones(mom_eq.n_elems)
		n = n_salt*to.ones(mom_eq.n_elems)
		# A = A_salt
		A = to.zeros(mom_eq.n_elems)
		A = A.double()
		T0_field = T0_field.double()
		ind_salt = to.tensor(ind_salt, dtype=to.long)
		A[ind_salt] = A_salt * to.exp(Q[0] / (R * T0_field[0]))
		A[ind_ovb] = 0.0
		creep_0 = sf.DislocationCreep(A, Q, n, "creep")
		A[ind_ovb] = 0.0
		creep_0 = sf.DislocationCreep(A, Q, n, "creep")

		# Create constitutive model
		mat.add_to_elastic(spring_0)
		# mat.add_to_non_elastic(kelvin)
		mat.add_to_non_elastic(creep_0)

		# Set constitutive model
		mom_eq.set_material(mat)

		# Set body forces
		g = -9.81
		g_vec = [0.0, 0.0, g]
		mom_eq.build_body_force(g_vec)

		# Time settings for equilibrium stage
		tc_eq = sf.TimeControllerParabolic(n_time_steps=20, initial_time=0.0, final_time=5, time_unit="day")
		# tc_eq = sf.TimeController(dt=0.5, final_time=5, initial_time=0.0, time_unit="hour")

		# Boundary conditions
		time_values = [0*ut.hour,  1*ut.hour]
		nt = len(time_values)

		bc_west_salt = momBC.DirichletBC(boundary_name="WEST", component=0, values=[0.0, 0.0], time_values=[0.0, tc_eq.t_final])
		bc_west_ovb = momBC.DirichletBC(boundary_name = "WEST_OB", component=0, values=[0.0, 0.0], time_values=[0.0, tc_eq.t_final])

		bc_east_salt = momBC.DirichletBC(boundary_name="EAST", component=0, values=[0.0, 0.0], time_values=[0.0, tc_eq.t_final])
		bc_east_ovb = momBC.DirichletBC(boundary_name = "EAST_OB", component=0, values=[0.0, 0.0], time_values=[0.0, tc_eq.t_final])

		bc_bottom = momBC.DirichletBC(boundary_name="BOTTOM", component=2, values=[0.0, 0.0], time_values=[0.0, tc_eq.t_final])

		bc_south_salt = momBC.DirichletBC(boundary_name="SOUTH", component=1, values=[0.0, 0.0], time_values=[0.0, tc_eq.t_final])
		bc_south_ovb = momBC.DirichletBC(boundary_name="SOUTH_OB", component=1, values=[0.0, 0.0], time_values=[0.0, tc_eq.t_final])

		bc_north_salt = momBC.DirichletBC(boundary_name="NORTH", component=1, values=[0.0, 0.0], time_values=[0.0, tc_eq.t_final])
		bc_north_ovb = momBC.DirichletBC(boundary_name="NORTH_OB", component=1, values=[0.0, 0.0], time_values=[0.0, tc_eq.t_final])

		# Extract geometry dimensions
		Lx = grid.Lx
		Ly = grid.Ly
		Lz = grid.Lz
		z_surface = 0.0

		g = 9.81
		ovb_thickness, salt_thickness, hanging_wall = get_geometry_parameters(grid_path)
		cavern_roof = ovb_thickness + hanging_wall
		p_roof = 0 + salt_density*g*hanging_wall + ovb_density*g*ovb_thickness

		# Pressure at the top of the salt layer (bottom of overburden)
		p_top = ovb_density*g*ovb_thickness

		bc_top = momBC.NeumannBC(boundary_name = "TOP",
							direction = 2,
							density = 0.0,
							ref_pos = z_surface,
							values = [0*MPa, 0*MPa],
							time_values = [0*day,  10*day],
							g = g_vec[2])

		bc_cavern = momBC.NeumannBC(boundary_name = "CAVERN",
							direction = 2,
							density = gas_density,
							ref_pos = cavern_roof,
							values = [P_max*p_roof, P_max*p_roof],
							time_values = [0*day,  10*day],
							g = g_vec[2])

		bc_equilibrium = momBC.BcHandler(mom_eq)
		bc_equilibrium.add_boundary_condition(bc_west_salt)
		bc_equilibrium.add_boundary_condition(bc_west_ovb)
		bc_equilibrium.add_boundary_condition(bc_east_salt)
		bc_equilibrium.add_boundary_condition(bc_east_ovb)
		bc_equilibrium.add_boundary_condition(bc_bottom)
		bc_equilibrium.add_boundary_condition(bc_south_salt)
		bc_equilibrium.add_boundary_condition(bc_south_ovb)
		bc_equilibrium.add_boundary_condition(bc_north_salt)
		bc_equilibrium.add_boundary_condition(bc_north_ovb)
		bc_equilibrium.add_boundary_condition(bc_top)
		bc_equilibrium.add_boundary_condition(bc_cavern)


		# Set boundary conditions
		mom_eq.set_boundary_conditions(bc_equilibrium)

		# Equilibrium output folder
		ouput_folder_equilibrium = os.path.join(output_folder, "equilibrium")

		# Print output folder
		if MPI.COMM_WORLD.rank == 0:
			print(ouput_folder_equilibrium)
			sys.stdout.flush()

		# Create output handlers
		output_mom = sf.SaveFields(mom_eq)
		output_mom.set_output_folder(ouput_folder_equilibrium)
		output_mom.add_output_field("u", "Displacement (m)")
		# output_mom.add_output_field("Temp", "Temperature (K)")
		output_mom.add_output_field("eps_tot", "Total strain (-)")
		output_mom.add_output_field("p_elems", "Mean stress (Pa)")
		output_mom.add_output_field("q_elems", "Von Mises stress (Pa)")
		output_mom.add_output_field("p_nodes", "Mean stress (Pa)")
		output_mom.add_output_field("q_nodes", "Von Mises stress (Pa)")
		output_mom.add_output_field("sig", "Stress (Pa)")
		outputs = [output_mom]

		# Define simulator
		sim = sf.Simulator_M(mom_eq, tc_eq, outputs, True)
		sim.run()

		# Print time
		if MPI.COMM_WORLD.rank == 0:
			end_time = MPI.Wtime()
			elaspsed_time = end_time - start_time
			formatted_time = time.strftime("%H:%M:%S", time.gmtime(elaspsed_time))
			print(f"Time: {formatted_time} ({elaspsed_time} seconds)\n")
			sys.stdout.flush()


		# # Time settings for operation stage
		# tc_op = sf.TimeController(dt=2, initial_time=0.0, final_time=240, time_unit="hour")
		# p_base = [0.8*p_roof, 0.2*p_roof, 0.2*p_roof, 0.8*p_roof, 0.8*p_roof]
		# t_values = [i*day for i in range(0, 10, 2)]
		# p_values = (p_base * (len(t_values) // len(p_base) + 1))[:len(t_values)]


		tc_op = sf.TimeController(dt=0.5, initial_time=0.0, final_time=100, time_unit="day")
		p_base = [P_max * p_roof, P_min * p_roof]
		t_values = [i * day for i in range(0, 100, 2)]
		p_values = (p_base * (len(t_values) // len(p_base) + 1))[:len(t_values)]

		# # Boundary conditions
		bc_west_salt = momBC.DirichletBC(boundary_name="WEST", component=0, values=[0.0, 0.0], time_values=[0.0, tc_op.t_final])
		bc_west_ovb = momBC.DirichletBC(boundary_name = "WEST_OB", component=0, values=[0.0, 0.0], time_values=[0.0, tc_op.t_final])

		bc_east_salt = momBC.DirichletBC(boundary_name="EAST", component=0, values=[0.0, 0.0], time_values=[0.0, tc_op.t_final])
		bc_east_ovb = momBC.DirichletBC(boundary_name = "EAST_OB", component=0, values=[0.0, 0.0], time_values=[0.0, tc_op.t_final])

		bc_bottom = momBC.DirichletBC(boundary_name="BOTTOM", component=2, values=[0.0, 0.0], time_values=[0.0, tc_op.t_final])

		bc_south_salt = momBC.DirichletBC(boundary_name="SOUTH", component=1, values=[0.0, 0.0], time_values=[0.0, tc_op.t_final])
		bc_south_ovb = momBC.DirichletBC(boundary_name="SOUTH_OB", component=1, values=[0.0, 0.0], time_values=[0.0, tc_op.t_final])

		bc_north_salt = momBC.DirichletBC(boundary_name="NORTH", component=1, values=[0.0, 0.0], time_values=[0.0, tc_op.t_final])
		bc_north_ovb = momBC.DirichletBC(boundary_name="NORTH_OB", component=1, values=[0.0, 0.0], time_values=[0.0, tc_op.t_final])

		bc_cavern = momBC.NeumannBC(
		    boundary_name="CAVERN",
		    direction=2,
		    density=gas_density,
		    ref_pos=cavern_roof,
		    values= p_values, 
		    time_values=t_values, 
		    g=g_vec[2]
		)


		bc_operation = momBC.BcHandler(mom_eq)
		bc_operation.add_boundary_condition(bc_west_salt)
		bc_operation.add_boundary_condition(bc_west_ovb)
		bc_operation.add_boundary_condition(bc_east_salt)
		bc_operation.add_boundary_condition(bc_east_ovb)
		bc_operation.add_boundary_condition(bc_bottom)
		bc_operation.add_boundary_condition(bc_south_salt)
		bc_operation.add_boundary_condition(bc_south_ovb)
		bc_operation.add_boundary_condition(bc_north_salt)
		bc_operation.add_boundary_condition(bc_north_ovb)
		bc_operation.add_boundary_condition(bc_top)
		bc_operation.add_boundary_condition(bc_cavern)

		# Set boundary conditions
		mom_eq.set_boundary_conditions(bc_operation)

		# Define output folder
		output_folder_operation = os.path.join(output_folder, "operation")

		# Print output folder
		if MPI.COMM_WORLD.rank == 0:
			print(output_folder_operation)
			sys.stdout.flush()

		# Create output handlers
		output_mom = sf.SaveFields(mom_eq)
		output_mom.set_output_folder(output_folder_operation)
		output_mom.add_output_field("u", "Displacement (m)")
		# output_mom.add_output_field("Temp", "Temperature (K)")
		output_mom.add_output_field("eps_tot", "Total strain (-)")
		output_mom.add_output_field("p_elems", "Mean stress (Pa)")
		output_mom.add_output_field("q_elems", "Von Mises stress (Pa)")
		output_mom.add_output_field("p_nodes", "Mean stress (Pa)")
		output_mom.add_output_field("q_nodes", "Von Mises stress (Pa)")
		outputs = [output_mom]

		# Define simulator
		sim = sf.Simulator_M(mom_eq, tc_op, outputs, False)
		sim.run()

		# # Print time
		# if MPI.COMM_WORLD.rank == 0:
		# 	end_time = MPI.Wtime()
		# 	elaspsed_time = end_time - start_time
		# 	formatted_time = time.strftime("%H:%M:%S", time.gmtime(elaspsed_time))
		# 	print(f"Time: {formatted_time} ({elaspsed_time} seconds)\n")
		# 	sys.stdout.flush()

		# save_cavern_pressure_to_csv(bc_cavern, output_folder)

		def save_cavern_pressure_to_csv(time_values, pressure_values, output_folder):
			output_dir = os.path.join(output_folder, "bc_cavern")
			os.makedirs(output_dir, exist_ok=True)
			file_path = os.path.join(output_dir, "cavern_pressure.csv")

			time_day = np.array(time_values) / day
			pressure_MPa = np.array(pressure_values) / MPa

			if MPI.COMM_WORLD.rank == 0:
				with open(file_path, "w") as f:
					f.write("time_day,pressure_MPa\n")
					for t, p in zip(time_day, pressure_MPa):
						f.write(f"{t},{p}\n")

				print(f"[info] Cavern pressure saved to: {file_path}")

		save_cavern_pressure_to_csv(t_values, p_values, output_folder_operation)


	main()	

