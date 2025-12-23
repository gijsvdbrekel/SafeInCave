import os
import sys
import time
from petsc4py import PETSc
sys.path.append(os.path.join("..", "..", "..", "safeincave"))
from build_input_file import create_input_file, MPa, GPa
from safeincave.Utils import read_json

# Start clock
start = time.time()

# Choose ouput folder name
results_folder = "case_2"

# Create input file
create_input_file(	
					P_min = 0.2, 
					P_max = 0.8,
					rho_ob = 2663,
					E_ob = 37*GPa,
					nu_ob = 0.19,
					rho_salt = 2219,
					E_salt = 87.21*GPa,
					nu_salt = 0.329,
					A_salt = 1.5e-38,
					n_salt = 4.3,
					case_folder = results_folder
)

# Stop clock
print(f"Total simulation time: {(time.time() - start)/60} minutes.")


