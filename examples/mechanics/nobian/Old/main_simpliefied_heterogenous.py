import os
import sys
import numpy as np
sys.path.append(os.path.join("..", "..", "safeincave"))
from InputFileAssistant import BuildInputFile
from Simulator import Simulator
from TimeSteppingDesigner import build_mapping
from TimeSteppingDesigner import wrap_cumulative
import matplotlib.pyplot as plt

# Useful units
hour = 60*60
day = 24*hour
year = 365*day
MPa = 1e6
GPa = 1e9

Equilibrium_toggle = True
PSC_toggle = False  # Toggle for pressure solution creep
P_soft_shutin = None #0.7  # Pressure threshold for soft shut-in (0.85 MPa)
cavern_fname = "_clay_grid"

# ---------- EQUILIBRIUM dynamic stepping -----------------------------
eq_min_dt  = 0.05*hour           # or what you like
eq_max_dt  = 4*hour	#1*day 	#3*day #5*day

t_max_eq = 	5*day	#30*day 	#100*day #200*day #0.2*day
eq_max_iterations = 500

eq_times   = [0*hour, 0.5*hour, t_max_eq]          # anchor pseudo-times 	100*day,
eq_dt      = [eq_min_dt, 0.5*hour,eq_max_dt] # target Δt values 		3*day, 


#HARD SHUTIN
# ---------- OPERATION dynamic stepping -------------------------------
min_time_step = 0.1*hour  # Minimum time step for the operation stage
max_time_step = 500*day #80*day  # Maximum time step for the operation stage

t_last = 300*year #20*year

times      = [ 0*day, 0.1*day, 1*day, 10*day, 50*day, 1*year, 5*year, 30*year, t_last]          # anchor times
desired_dt = [min_time_step, 0.5*hour, 2*hour, 15*hour, 20*hour, 10*day, 60*day, 150*day, max_time_step]       # target step sizes


# #SOFT SHUTIN
# # ---------- OPERATION dynamic stepping -------------------------------
# max_time_step = 8*day  # Maximum time step for the operation stage
# min_time_step = 0.1*hour  # Minimum time step for the operation stage
# t_last = 20*year

# # CONFIGURATION FOR TIME MAPPING with anchor points
# times      = [ 0*day, 0.1*day, 1*day, 10*day, 125*day, 365*day, t_last]          # anchor times
# desired_dt = [min_time_step, 0.5*hour, 1*hour, 8*hour, 4*day, 6*day, max_time_step]       # target step sizes


# # -CONFIGURATION FOR MAPPING WITH CUMULATIVE FUNCTION
# tau   = 6*hour
# f_cum = lambda t: t + 0.5*(t**2)/tau          # quadratic tail
# fun_fun = wrap_cumulative(f_cum,
#                                t_final = t_last,
#                                scale   = min_time_step,   # 1800 s
#                                h       = 1.0)              # FD stencil (s)

def get_geometry_parameters(path_to_grid):
	f = open(os.path.join(path_to_grid, "geom.geo"), "r")
	data = f.readlines()
	ovb_thickness = float(data[10][len("ovb_thickness = "):-2])
	salt_thickness = float(data[11][len("salt_thickness = "):-2])
	hanging_wall = float(data[12][len("hanging_wall = "):-2])
	return ovb_thickness, salt_thickness, hanging_wall

def create_input_file():
	# Initialize input file object
	ifa = BuildInputFile()

	# Create input_grid section
	path_to_grid = os.path.join("..", "..", "grids_final", cavern_fname)
	ifa.set_input_grid(path_to_grid, "geom")

	# Create output section
	ifa.set_output_folder(os.path.join("output", "case_0"))

	# Create solver settings section
	# ifa.set_krylov_solver(method="cg", preconditioner="petsc_amg", rel_tol=1e-12)
	ifa.set_krylov_solver(method="cg", preconditioner="ilu", rel_tol=1e-12)
	# ifa.set_krylov_solver(method="bicgstab", preconditioner="ilu", rel_tol=1e-12)
	# ifa.set_krylov_solver(method="bicgstab", preconditioner="petsc_amg", rel_tol=1e-12)
	# ifa.set_direct_solver(method="petsc")

	# Create simulation_settings section
	ifa.set_equilibrium_stage(active=Equilibrium_toggle, dt=0.02*day, tol=1e-9, ite_max=eq_max_iterations, max_time_step=eq_max_dt, min_time_step=eq_min_dt, t_max=t_max_eq, n_skip=3)
	ifa.set_operation_stage(active=True, dt=3*hour, n_skip=10, max_time_step=max_time_step, min_time_step=min_time_step)

	# Define densities
	salt_density = 2200
	salt_top_density = salt_density  # Density at the top of the salt layer
	salt_bottom_density = salt_density  # Density at the bottom of the salt layer
	ovb_density = 2800
	clay_density = 2500  # Density of the clay layer
	
	gas_density = 10


	# Create body_forces section
	#Order in .geo file is: overburden, salt_top, interlayer, salt_bottom
	##IT USED TO BE SALT, OVB
	ifa.section_body_forces(density=[ovb_density, salt_top_density, clay_density, salt_bottom_density], direction=2)

	# Create time_settings section
	time_list = [0*day, t_last]
	ifa.section_time(time_list, theta=0.5)

	# Create boundary_conditions section

	# Add Dirichlet boundary conditions
	ifa.add_dirichlet(name="West", values=list(np.zeros(len(time_list))), component=0)
	ifa.add_dirichlet(name="South", values=list(np.zeros(len(time_list))), component=1)
	ifa.add_dirichlet(name="Bottom", values=list(np.zeros(len(time_list))), component=2)
	ifa.add_dirichlet(name="East", values=list(np.zeros(len(time_list))), component=0)
	ifa.add_dirichlet(name="North", values=list(np.zeros(len(time_list))), component=1)

	# Extract geometry dimensions
	Lx = ifa.grid.Lx
	Ly = ifa.grid.Ly
	Lz = ifa.grid.Lz
	z_surface = 0.0

	g = 9.81
	ovb_thickness, salt_thickness, hanging_wall = get_geometry_parameters(path_to_grid)
	cavern_roof = ovb_thickness + hanging_wall
	p_roof = 0 + salt_density*g*hanging_wall + ovb_density*g*ovb_thickness
	# Pressure at the top of the salt layer (bottom of overburden)
	p_top = ovb_density*g*ovb_thickness


	#pseudocode for how simulation name to be set
	# Format physical values into readable suffixes
	depth_m = int(round(cavern_roof))  # Round to integer meters
	duration_y = int(round(t_last / year))  # Convert seconds to years
	eq_days = int(round(t_max_eq / day))  # Convert seconds to days
	soft_shut_str = f"{P_soft_shutin*100:.2f}" if P_soft_shutin is not None else "no"
	psc_str = "on" if PSC_toggle else "off"

	simulation_name = f"{depth_m}m_{duration_y}y_{eq_days}d_soft-{soft_shut_str}_PSC-{psc_str}"

	ifa.pass_cavern_parameters(depth=cavern_roof, P_lithosphere=p_roof, simulation_name=simulation_name)

	# Add Neumann boundary condition
	p_wellhead= 1168.5291401077016 * 9.81 * cavern_roof  # Pressure at the wellhead (in Pa)
	print(p_roof, cavern_roof, p_wellhead)

	if P_soft_shutin is not None:
		P_soft_shutin_real = P_soft_shutin * p_roof
		ifa.add_special_neumann(name="Cavern", P_0=p_wellhead, T_0=273+40, direction=2, density=1168.5291401077016, reference_position=cavern_roof, P_threshold=P_soft_shutin_real) #, P_threshold=0.85*p_roof
	else:
		ifa.add_special_neumann(name="Cavern", P_0=p_wellhead, T_0=273+40, direction=2, density=1168.5291401077016, reference_position=cavern_roof, P_threshold=None)
	# ifa.add_neumann(name="Cavern", values=[0.8*p_roof, 0.8*p_roof, 0.8*p_roof, 0.8*p_roof, 0.8*p_roof], direction=2, density=1000, reference_position=cavern_roof)
	ifa.add_neumann(name="Top", values=[0*MPa]*len(time_list), direction=2, density=0.0, reference_position=z_surface)
	# ifa.add_special_neumann(name="Top", P_0=0*MPa, T_0=273+40, direction=2, density=0.0, reference_position=z_surface)


	# Define constitutive model
	# Add elastic elements
	ifa.add_elastic_element(
							name="Spring_0",
							E=[15*GPa, 102*GPa, 8*GPa, 102*GPa],  # E for overburden, salt top, interlayer, salt bottom
							nu=0.32,
							active=True,
							equilibrium=True
	)

	# Add viscoelastic elements
	ifa.add_viscoelastic_element(
							name="KelvinVoigt_0",
							E=10*GPa,
							nu=0.32,
							eta=105e11,
							active=False,
							equilibrium=False
	)

	# Build geothermal profile
	def geothermal_grad(x, y, z):
		km = 1000
		dTdZ = 27/km
		T_surface = 20 + 273
		return T_surface - dTdZ*z
	T_profile = ifa.build_custom_field(geothermal_grad)

	# Add inelastic elements
	ifa.add_dislocation_creep_element(
							name="disCreep",
							A=[0.0, 1.1e-21, 0, 1.1e-21], # A for overburden, salt top, interlayer, salt bottom
							n=3,#n=3.0,
							Q=51600, #Q=51600,
							T=T_profile,
							active=True,
							equilibrium=True
	)

	def Aps(x, y, z):
		if z < -400:
			km = 1000
			dTdZ = 27/km
			T_surface = 20 + 273
			T = T_surface - dTdZ*z
			B_ps = 1.29e-13
			A_ps = B_ps/(273 + T)
		else:
			A_ps = 0.0
		return A_ps

	A_ps = ifa.build_custom_field(Aps)
	ifa.add_dislocation_creep_element(
							name="PressureSolution",
							A=A_ps,
							n=1.0,
							Q=13184,
							T=T_profile,
							active=PSC_toggle,
							equilibrium=PSC_toggle
	)

	ifa.add_desai_element(
							name="desai",
							mu_1=[5.3665857009859815e-11, 0.0],
							N_1=3.1,
							n=3.0,
							a_1=1.965018496922832e-05,
							eta=0.8275682807874163,
							beta_1=0.0048,
							beta=0.995,
							m=-0.5,
							gamma=0.095,
							alpha_0=0.0022,
							sigma_t=5.0,
							active=False,
							equilibrium=False
	)


	# Save input_file.json
	ifa.save_input_file("input_file.json")

	return ifa.input_file


def main():
	# Create input file
	input_file = create_input_file()
	
	# Build simulator
	sim = Simulator(input_file)

	op_time_map = build_mapping(times, desired_dt, method="pchip")
	eq_time_map = build_mapping(eq_times, eq_dt, method="pchip")

	# n_plot   = 400
	# t_plot   = np.linspace(0.0, t_last, n_plot)
	# g_plot   = op_time_map(t_plot)


	# # Plot the time mapping NOTE: MUST CLOSE PLOT WINDOW TO CONTINUE
	# plt.figure()
	# plt.plot(t_plot / 3600.0, g_plot / 3600.0, '.-')
	# plt.xlabel("Pseudo-time t  [h]")
	# plt.ylabel("Step size g(t)  [h]")
	# plt.title("Dynamic time-step profile")
	# plt.grid(True, ls=":")
	# plt.tight_layout()
	# plt.show()

	# time_map = fun_fun  # Use the custom time mapping function
	sim.set_equilibrium_time_transform_function(eq_time_map)
	sim.set_time_transform_function(op_time_map)

	# Run simulation
	sim.run()

if __name__ == '__main__':
	main()



