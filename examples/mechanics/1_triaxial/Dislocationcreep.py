import numpy as np

R = 8.314  # J/mol/K
T = 298.0  # K

# stress range in Pa (1..100 MPa)
sigma = np.logspace(6, 8, 200)

# OLD
A_old = 40.0 * (1e-6)**4.6 / (365.25*24*3600)   # 1/s / (Pa^n)  (?) check!
Q_old = 6495.0 * 8.32                            # J/mol
n_old = 4.6

# NEW
A_new = 1.9e-20                                  # check units!
Q_new = 51600.0                                  # J/mol
n_new = 3.0

epsdot_old = A_old * (sigma**n_old) * np.exp(-Q_old/(R*T))
epsdot_new = A_new * (sigma**n_new) * np.exp(-Q_new/(R*T))

x = np.log10(sigma/1e6)        # log10(MPa)
y_old = np.log10(epsdot_old)
y_new = np.log10(epsdot_new)

print("Old @10MPa:", epsdot_old[np.argmin(np.abs(sigma-10e6))])
print("New @10MPa:", epsdot_new[np.argmin(np.abs(sigma-10e6))])
