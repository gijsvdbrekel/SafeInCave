import itasca as it
import numpy as np
import pylab as plt
from scipy.interpolate import interp1d
import pandas as pd

e_axial = None # axial strain
S_diff = None # different stress (s3 - s1)
e_vol = None # dilatancy
e_rate = None # axial strain rate
time = None # time as a function of step cycle 

def get_data_compression():
    i = it.cycle()
    # collect every 10 steps
    if i %10 != 0:
        return
    zn = it.zone.find(1)
    gp = it.gridpoint.find(8)
    e_axial.append(-gp.disp_z()*100)
    S_diff.append((zn.stress()[0, 0] - zn.stress()[2, 2])*1e-6)
    e_vol.append(zn.strain_vol_inc()*100)

    
def get_data_extension():
    i = it.cycle()
    # collect every 10 steps
    if i %10 != 0:
        return
    zn = it.zone.find(1)
    gp = it.gridpoint.find(8)
    e_axial.append(-gp.disp_z()*100)
    S_diff.append((zn.stress()[0, 0] - zn.stress()[2, 2])*1e-6)
    e_vol.append(zn.strain_vol_inc()*100)
    e_rate.append(zn.strain_rate()[2, 2])
    time.append(i)    


def update_stress():
    """
    Deprecated method
    """
    i = it.cycle()
    global GSa
    global GdSr
    val = GSa + GdSr*i
    it.command("zone face apply-remove stress-normal range group 'North' or 'South' or 'East' or 'West' ")
    it.command(f"zone face apply stress-normal {val} range group 'North' or 'South' or 'East' or 'West' ")


def compression_model(variables:np.array, model:str, Sr:float, eps:float, nstep:int, e_max:float):
    """
    Run a simulation of the digital twin of a compression triaxial test with a define mechanical  model.
    
    Parameters:
    - variables: numpy vector, dependent of the model chosen:
        - mohr-coulomb: 
            - 0: young modulus [Pa]
            - 1: poisson ratio
            - 2:cohesion [Pa]
            - 3:friction angle [°]
            - 4: dilation [°])
        - strain-soft: see code
        - plastic-hardening: see code
    - model: str : mechanical model used
    - Sr: initial Confining stress in the triaxial test [Pa]
    - eps: fixed Strain rate during the experiment [s-1]
    - nstep: number of step. In the definition of the model is correspond to the model duration [s] as the strain rate correspond to one step.
    - e_max: maximum axial strain needed. This will define for how much steps the model is run
    
    Return: return a hisotry of the system deformation
    - e_axial: history of axial strain [%]
    - S_diff: differential stress [Pa]
    - e_vol: volumetrix strain or dilatancy [%]
    
    """
    
    it.command(f"""
    model new
    model large-strain off
    
    zone create brick size 1 1 1
    zone face skin
    
    zone cmodel {model}
    """)

    global e_axial
    e_axial = []
    global S_diff
    S_diff = []
    global e_vol
    e_vol = []
    global e_rate
    
    if model == "mohr-coulomb" and variables.shape[0] == 5:
        it.command(f"zone property young={variables[0]} poisson={variables[1]} cohesion={variables[2]} friction={variables[3]} dilation={variables[4]}")
    elif model == "strain-soft" and variables.shape[0] == 10:
        it.command(f"""
        zone property young={variables[0]} poisson={variables[1]} cohesion={variables[2]} friction={variables[3]}
        zone property table-friction '1'
        zone property table-dilation '2'
        table '1' add (0.0,variables[4]), (0.05,variables[5]),(0.10,variables[6])
        table '2' add (0.0,variables[7]), (0.05,variables[8]),(0.10,variables[9])
        """)
    elif model == "plastic-hardening" and variables.shape[0] == 10:
        it.command(f"""
        zone property stiffness-50-reference={variables[0]} stiffness-ur-reference={variables[1]} pressure-reference={variables[2]} exponent={variables[3]}
        zone property cohesion={variables[4]} friction={variables[5]} dilation={variables[6]}
        zone property stress-1-effective={variables[7]} stress-2-effective={variables[8]} stress-3-effective={variables[9]}
        
        """)
    
    else:
        raise Exception(f"Model {model} with {varaibles.shape[0]} variables is not supported")
    
    
    it.command(f"""
    zone initialize stress xx {Sr} yy {Sr} zz {Sr}
    zone face apply velocity-z 0.0 range group "Bottom"
    zone face apply velocity-z {eps} range group "Top"
    zone face apply stress-normal {Sr} range group "North" or "South" or "East" or "West"
    """)
    
    
    it.set_callback("get_data_compression", 11.0)
    
    it.command(f"model step {nstep}")
    
    while e_axial[-1] < e_max:
        it.command(f"model step {nstep}")
    
    return e_axial, S_diff, e_vol
    
    
def extension_model(variables:np.array, model:str, Sa:float, dSr:float, holdstep:int, nstep:int):
    """
    Extension test model in which the confining pressure is solvely increase while the axial pressure remains the same
    
    Parameters:
    - variables: variable for the material model. see compression test.
    - model: type of model use. Supported are mohr-coulomb, plastic-hardening, strain-softening
    - Sa: initial confining pressure and axial stress in Pa
    - dSr: rate of radial stress increase. 
    - holdstep: number of steps where the stress is held fixed.
    - nstep: number of step where the stress is build up.
    
    """
    it.command(f"""
    model new
    model large-strain off
    
    zone create brick size 1 1 1
    zone face skin
    
    zone cmodel {model}
    """)
    
    it.fish.set('GSa', Sa)
    it.fish.set('GdSr', dSr)
    global e_axial
    e_axial = []
    global S_diff
    S_diff = []
    global e_vol
    e_vol = []    
    global e_rate
    e_rate = []
    global time
    time = []
    
    if model == "mohr-coulomb" and variables.shape[0] == 5:
            
        it.command(f"zone property young={variables[0]} poisson={variables[1]} cohesion={variables[2]} friction={variables[3]} dilation={variables[4]}")
    elif model == "strain-soft" and variables.shape[0] == 10:
        it.command(f"""
        zone property young={variables[0]} poisson={variables[1]} cohesion={variables[2]} friction={variables[3]}
        zone property table-friction '1'
        zone property table-dilation '2'
        table '1' add (0.0,variables[4]), (0.05,variables[5]),(0.10,variables[6])
        table '2' add (0.0,variables[7]), (0.05,variables[8]),(0.10,variables[9])
        """)
    elif model == "plastic-hardening" and variables.shape[0] == 10:
        it.command(f"""
        zone property stiffness-50-reference={variables[0]} stiffness-ur-reference={variables[1]} pressure-reference={variables[2]} exponent={variables[3]}
        zone property cohesion={variables[4]} friction={variables[5]} dilation={variables[6]}
        zone property stress-1-effective={variables[7]} stress-2-effective={variables[8]} stress-3-effective={variables[9]}
        
        """)
    
    else:
        raise Exception(f"Model {model} with {varaibles.shape[0]} variables is not supported")
    
    it.command(f"""
    fish define stress_bc(zn, i)
        k = mech.cycle()
        if k<= {holdstep}
            return 1.
        endif
        
        global GSa
        global GdSr
        val = 1 + GdSr*(k - {holdstep})/GSa
        ;io.out('Stress: ' + string(val) + ' cycles ' + string(k))
        return val
    end
    
    """)
    
    it.command(f"""
    zone initialize stress xx {Sa} yy {Sa} zz {Sa}
    zone face apply velocity-z 0.0 range group "Bottom"
    zone face apply stress-normal {Sa} range group "Top"
    """)

    it.command(f"zone face apply stress-normal {Sa} fish-local stress_bc range group 'North' or 'South' or 'East' or 'West' ")

    it.set_callback("get_data_extension", 11.0)

    it.command(f"model step {holdstep+ nstep}")
    
    return None    
    
def extension_model_old(variables:np.array, model:str, Sa:float, dSr:float, nstep:int):

    it.command(f"""
    model new
    model large-strain off
    
    zone create brick size 1 1 1
    zone face skin
    
    zone cmodel {model}
    """)
    global GSa
    GSa = Sa
    global GdSr
    GSr = dSr
    global e_axial
    e_axial = []
    global S_diff
    S_diff = []
    global e_vol
    e_vol = []    
    
    if model == "mohr-coulomb" and variables.shape[0] == 5:
            
        it.command(f"zone property young={variables[0]} poisson={variables[1]} cohesion={variables[2]} friction={variables[3]} dilation={variables[4]}")
    elif model == "strain-soft" and variables.shape[0] == 10:
        it.command(f"""
        zone property young={variables[0]} poisson={variables[1]} cohesion={variables[2]} friction={variables[3]}
        zone property table-friction '1'
        zone property table-dilation '2'
        table '1' add (0.0,variables[4]), (0.05,variables[5]),(0.10,variables[6])
        table '2' add (0.0,variables[7]), (0.05,variables[8]),(0.10,variables[9])
        """)
    elif model == "plastic-hardening" and variables.shape[0] == 10:
        it.command(f"""
        zone property stiffness-50-reference={variables[0]} stiffness-ur-reference={variables[1]} pressure-reference={variables[2]} exponent={variables[3]}
        zone property cohesion={variables[4]} friction={variables[5]} dilation={variables[6]}
        zone property stress-1-effective={variables[7]} stress-2-effective={variables[8]} stress-3-effective={variables[9]}
        
        """)
    
    else:
        raise Exception(f"Model {model} with {varaibles.shape[0]} variables is not supported")
    
    
    it.command(f"""
    zone initialize stress xx {Sa} yy {Sa} zz {Sa}
    zone face apply velocity-z 0.0 range group "Bottom"
    zone face apply stress-normal {Sa} range group "Top"
    """)
    
    it.set_callback("get_data", 11.0)
    it.set_callback("update_stress", -100.0)
    it.command(f"model step {nstep}")
    
    return None

def init_minimization_from_IfG_data(filename, config):
    """
    Extract configuration data from excel file provided by IfG
    This method is specific to IfG data
    """
    for k, v in config.items():
        df = pd.read_excel(filename, skiprows=[1], sheet_name=k)
        v["e_axial"] = df["Eps_axial"].to_numpy()
        v["s_diff"] = df["stress_diff"].to_numpy()
        v["e_vol"] = df["Eps_vol"].to_numpy()
        
    
    return config
    
    
def compression_operand(variables:np.array, model:str, data:dict, N:int):
    """
    This is a single operand used to be used by the optimization.
    It work with the data of a triaxial compression test experiment.
    
    Parameter:
    - variable to run the system on
    - data: dictionary of experiment parameters 
    - N: N number of samples to the curve for comparison
   
    """
    e_min = 0.
    e_max = data["e_max"]
    x_sample = np.linspace(np.sqrt(e_min + (e_max-e_min)/(N+1)), np.sqrt(e_max - (e_max-e_min)/(N+1)), N)**2
    
    
    m_e_axial, m_s_diff, m_e_vol = compression_model(variables, model, -data["s3"], -data["e_rate"], 2000, e_max)
    
    m_fit_stress = np.interp(x_sample, m_e_axial, m_s_diff, left=0.)
    m_fit_dilatancy = np.interp(x_sample, m_e_axial, m_e_vol, left=0.)
    d_fit_stress = np.interp(x_sample, data["e_axial"], data["s_diff"], left=0.)
    d_fit_dilatancy = np.interp(x_sample, data["e_axial"], data["e_vol"], left=0.)
    
    res = np.hstack([(d_fit_stress - m_fit_stress)/d_fit_stress, (d_fit_dilatancy - m_fit_dilatancy)/d_fit_dilatancy])
    
    return res

    
    
    
### actual code to run
#x0 = np.array([0.5e9, 0.3, 0.3e6, 30, 20])
#x0 = np.array([0.2e9, 0.3, 0.8e6, 30, 20, 25, 22, 0, 20, 1])
#x0 = np.array([0.2e8, 8.0e8, 1e5, 0.8, 0.3e6, 30, 10, -1e5, -1e5, -1e5])
#compression_model(x0, "mohr-coulomb", -1.5e6, -1e-5, 5000)
#extension_model(x0, "mohr-coulomb", -1.1e6, -6e6/3600, 2000)
