#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import geopandas as gpd
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def plot_SIRWellMixed(params,outputs):
    """
    plot SIRWellMixed model results
    S,I, R populations over time
    """	
    nbTimeSteps = len(outputs)
    time=np.zeros(nbTimeSteps)
    x = np.zeros([nbTimeSteps, 3])
    for i in range(0,nbTimeSteps):
        time[i] = outputs[i][0]
        x[i,:] = outputs[i][1]	
        
    plt.plot(time,x[:,0], 'red', label='Susceptible')
    plt.plot(time,x[:,1], 'green', label='Infected')
    plt.plot(time,x[:,2], 'blue', label='Recovered')
    plt.legend()
    string_for_title='SIRWellMixed, k_SI ='+str(params['k_SI']) +', k_IR = '+ str(params['k_IR']) + ', S_ini = '+ str(params['S_ini']) + ', I_ini = '+str( params['I_ini']) + ', R_ini = '+ str(params['R_ini'])
    plt.title(string_for_title)
    string_for_savefig = 'SIRWellMixed_k_SI='+str(params['k_SI']) +'_k_IR='+ str(params['k_IR']) + '_S_ini='+ str(params['S_ini']) + '_I_ini='+str( params['I_ini']) + '_R_ini='+ str(params['R_ini'])
    plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
    # plt.show()
	
def plot_SIRonRegions(params,outputs):
	"""
	TODO
	"""	

def plot_SIRonDepartments(params,outputs):
	"""
	TODO
	"""	
	
