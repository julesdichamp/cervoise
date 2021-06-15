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

def modelODEs_SIRLD(z, t, params):
    # model from Marc Lavielle: http://webpopix.org/covidix19.html
    # order of states x = S,  I, R, L, D
    k_SI_rescaled = params['k_SI'] / params['total_population']
    dSdt = - k_SI_rescaled * z[0] * z[1]
    dIdt =   k_SI_rescaled * z[0] * z[1] - params['k_IR'] * z[1] - params['k_IL'] * z[1]
    dRdt =   params['k_IR'] * z[1] 
    dLdt =   params['k_IL'] * z[1] - params['k_LD'] * z[3]
    dDdt =   params['k_LD'] * z[3]
    dxdt = [dSdt, dIdt, dRdt, dLdt, dDdt]
    return dxdt

def modelODEs_SIRLDmodified(z, t, params):
    # model from Marc Lavielle: http://webpopix.org/covidix19.html
    # order of states x = S,  I, R, L, D
    k_SI_rescaled = params['k_SI'] / params['total_population']
    dSdt = - k_SI_rescaled * z[0] * z[1]
    dIdt =   k_SI_rescaled * z[0] * z[1] - params['k_IR'] * z[1] - params['k_IL'] * z[1]
    dRdt =   params['k_IR'] * z[1] + params['k_LR'] * z[3]  
    dLdt =   params['k_IL'] * z[1] - params['k_LR'] * z[3] - params['k_LD'] * z[3]
    dDdt =   params['k_LD'] * z[3]
    dxdt = [dSdt, dIdt, dRdt, dLdt, dDdt]
    return dxdt

def modelODEs_SIR(z, t, params):
    # order of states x = S, I, R
    k_SI_rescaled = params['k_SI'] / params['S_ini']
    dSdt = - k_SI_rescaled * z[0] * z[1]
    dIdt =   k_SI_rescaled * z[0] * z[1] - params['k_IR'] * z[1] 
    dRdt =   params['k_IR'] * z[1]
    dxdt = [dSdt, dIdt, dRdt]
    return dxdt

def modelODEs_SEIsIasHoRadRD(z, t, params):
    # order of states x = S, E, Is, I_as, Ho, R, Rad, D
    k_SE_rescaled = params['k_SE'] / params['total_population']
    dSdt   = - k_SE_rescaled * z[0] * (z[2] + z[3])
    dEdt   =   k_SE_rescaled * z[0] * (z[2] + z[3]) - params['k_EIas'] * z[1] - params['k_EIs'] * z[1]
    dIasdt =   params['k_EIas'] * z[1] - params['k_IasR'] * z[2]
    dIsdt  =   params['k_EIs'] * z[1]  - params['k_IsR'] * z[3] - params['k_IsHo'] * z[3]
    dHodt  =   params['k_IsHo'] * z[3] - params['k_HoD'] * z[4] - params['k_HoRad'] * z[4]
    dRdt   =   params['k_IsR'] * z[3] + params['k_IasR'] * z[2]
    dRaddt =   params['k_HoRad'] * z[4]
    dDdt   =   params['k_HoD'] * z[4]
    dxdt   = [dSdt, dEdt, dIasdt, dIsdt, dHodt, dRdt, dRaddt, dDdt]
    return dxdt

def InitModelODEs(params, G):    
    params, G = InitModels(params, G)     
    time_ODEs =  np.linspace(params['time_start'], params['time_end'], 10000)    
    params['time_ODEs'] = time_ODEs
        
    return params , G
            
def LaunchModelODEs(params):
    
    if (params["model_name"] == "RealSEIRLDWellMixed"):
        x = odeint(modelODEs_SEIRLD, params['x_ini'], params['time_ODEs'], args=(params,))
    elif (params["model_name"] == "SIRWellMixed"):
        x = odeint(modelODEs_SIR, params['x_ini'], params['time_ODEs'], args=(params,))
    elif (params["model_name"] == "RealSIRLDWellMixed"):
        x = odeint(modelODEs_SIRLD, params['x_ini'], params['time_ODEs'], args=(params,))
    elif (params["model_name"] == "RealSIRLDmodifiedWellMixed"):
        x = odeint(modelODEs_SIRLDmodified, params['x_ini'], params['time_ODEs'], args=(params,))
    elif (params["model_name"] == "SIRWellMixed"):
        x = odeint(modelODEs_SIR, params['x_ini'], params['time_ODEs'], args=(params,))
    elif (params["model_name"] == "RealSEIsIasHoRadRDWellMixed"):
        x = odeint(modelODEs_SEIsIasHoRadRD, params['x_ini'], params['time_ODEs'], args=(params,))
        
    # write as outputs from Master equation framework in models.py
    outputs = []
    for i in range(0,len(params['time_ODEs'])):
        outputs.append([params['time_ODEs'][i], x[i,:]])
    
    return outputs
