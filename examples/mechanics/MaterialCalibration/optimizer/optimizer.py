from models import *
from scipy.optimize import minimize
import pylab as plt

def jacobian_1(x0:np.array, data:dict, model:str):
    """
    Calculate the jacobian of minimization_function_1 for the given array of variables
    
    Parameters
    - x0: variable values at which the jacobian is calculated
    - data: data for the model
    - model: type of model used
    """
    pass
    

def minimization_function_1(init:np.array, bounds, data:dict, model:str, maxiter):
    """
    This minimize the function based on the compression data loaded from a file. 
    
    A custom minimization function is created and a lot of values are hardcoded
    
    Parameters
    - init: initial guess for parameters
    - bounds: bounds of optimization
    - conf: data and model configuration to fit
    - model: model to optimise in flac
    - maxiter: maximum number of iterations
    """
    
    def objective_fct(x):
        
        objectives = []
        for k, v in data.items():
            objectives.append(np.sqrt(np.sum(compression_operand(x, model, v, 20)**2)))
        res = sum(objectives)
        print(f"--------------\n res: {res} variables: {x}")
        print(f"objective: \n {objectives} \n")
        return res
        
    result = minimize(objective_fct, init, method="L-BFGS-B", bounds=bounds, options={"maxiter":maxiter, "disp":True}) 
    
    return result

def minimization_function_2(init:np.array, bounds, data:dict, model:str, maxiter):
    """
    This minimize the function based on the compression data loaded from a file. 
    
    A custom minimization function is created and a lot of values are hardcoded
    
    Parameters
    - init: initial guess for parameters
    - bounds: bounds of optimization
    - conf: data and model configuration to fit
    - model: model to optimise in flac
    - maxiter: maximum number of iterations
    """
    
    def objective_fct(x):
        
        objectives = []
        try:
            for k, v in data.items():
                objectives.append(np.sqrt(np.sum(compression_operand(x, model, v, 20)**2)))
        except ValueError:
            res = np.nan
        else:
            res = sum(objectives)
        print(f"--------------\n res: {res} variables: {x}")
        print(f"objective: \n {objectives} \n")
        return res
        
    result = minimize(objective_fct, init, method="Nelder-Mead", bounds=bounds, options={"maxiter":maxiter, "disp":True}) 
    
    return result

def plot_stress_diff(x0:np.array, data:dict, model:str):
    """
    Plot the model against all data for the stress_difference
    """
    fig, axs = plt.subplots(2)
    e_max = 0.
    for k, v in data.items():
        axs[0].plot(v["e_axial"], v["s_diff"], label=k)
        axs[1].plot(v["e_axial"], v["e_vol"], label=k)
        e_axial, s_diff, e_vol = compression_model(x0, model, -v["s3"], -v["e_rate"], 2000, v["e_max"])
        axs[0].plot(e_axial, s_diff, label=f"model P={v['s3']*1e-6} MPa - e={v['e_rate']}", linestyle='dashed')
        axs[1].plot(e_axial, e_vol, label=f"model P={v['s3']*1e-6} MPa - e={v['e_rate']}", linestyle='dashed')
    axs[0].legend(loc="upper right")
    axs[1].legend(loc="upper right")
    axs[0].set_ylabel("Stress difference [MPa]")
    axs[1].set_ylabel("Dilatancy [%]")
    axs[1].set_xlabel("axial strain [%]")
    