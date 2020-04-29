#!/usr/bin/env python3

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from models import InitModels

def modelODEs_SEIRLD(z, t, params):    
    # order of states x = S, E, I, R, L, D
    k_SE_rescaled = params['k_SE'] / params['total_population']
    dSdt = - k_SE_rescaled * z[0] * z[2]
    dEdt =   k_SE_rescaled * z[0] * z[2] - params['k_EI'] * z[1]
    dIdt =   params['k_EI'] * z[1] - params['k_IR'] * z[2] - params['k_IL'] * z[2]
    dRdt =   params['k_IR'] * z[2] + params['k_LR'] * z[4]
    dLdt =   params['k_IL'] * z[2] - params['k_LR'] * z[4] - params['k_LD'] * z[4]
    dDdt =   params['k_LD'] * z[4]
    dxdt = [dSdt, dEdt, dIdt, dRdt, dLdt, dDdt]
    return dxdt

def InitModelODEs(params, G):
    if (params["model_name"] == "RealSEIRLDWellMixed"):        
        params, G = InitModels(params, G)     
        time_ODEs =  np.linspace(params['time_start'], params['time_end'], 10000)    
        params['time_ODEs'] = time_ODEs
        
    return params , G
            
def LaunchModelODEs(params):
    if (params["model_name"] == "RealSEIRLDWellMixed"):
        x = odeint(modelODEs_SEIRLD, params['x_ini'], params['time_ODEs'], args=(params,))
    
    # write as outputs from Master equation framework in models.py
    outputs = []
    for i in range(0,len(params['time_ODEs'])):
        outputs.append([params['time_ODEs'][i], x[i,:]])
    
    return outputs
