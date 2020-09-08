#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import geopandas as gpd
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from itertools import count
import copy
import seaborn as sns
import time
sns.set_context("talk")

def plot_SIRWellMixed(params,outputs):    
    """
    plot SIRWellMixed model results
    S,I, R populations over time
    """
    # save outputs
    f = open(params['outputs_dir']+'/outputs.txt',"w")
    f.write(str(outputs))
    f.close()
    
    # save params inputs
    f = open(params['outputs_dir']+'/params.txt',"w")
    f.write(str(params))
    f.close()
    
    nbTimeSteps = len(outputs)
    time=np.zeros(nbTimeSteps)
    x = np.zeros([nbTimeSteps, 3])
    for i in range(0,nbTimeSteps):
        time[i] = outputs[i][0]
        x[i,:] = outputs[i][1]

    fig2 = plt.figure(figsize=(16, 10)) # in inches!
    plt.plot(time,x[:,0], 'red', label='Susceptible')
    plt.plot(time,x[:,1], 'green', label='Infected')
    plt.plot(time,x[:,2], 'blue', label='Recovered')
    plt.legend()
    plt.xlabel('Arbitrary time unit')
    plt.ylabel('Population numbers')
    string_for_title='SIRWellMixed, k_SI ='+str(params['k_SI']) +', k_IR = '+ str(params['k_IR']) + ', S_ini = '+ str(params['S_ini']) + ', I_ini = '+str( params['I_ini']) + ', R_ini = '+ str(params['R_ini'])
    plt.title(string_for_title)
    string_for_savefig = 'SIRWellMixed_k_SI='+str(params['k_SI']) +'_k_IR='+ str(params['k_IR']) + '_S_ini='+ str(params['S_ini']) + '_I_ini='+str( params['I_ini']) + '_R_ini='+ str(params['R_ini'])
    plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
    # plt.close()
    plt.show()
def plot_SIRHigherOrderWellMixed(params,outputs):
    """
    plot SIRHigherOrderWellMixed model results
    S,I, R populations over time
    """
    # save outputs
    f = open(params['outputs_dir']+'/outputs.txt',"w")
    f.write(str(outputs))
    f.close()
    
    # save params inputs
    f = open(params['outputs_dir']+'/params.txt',"w")
    f.write(str(params))
    f.close()
    
    nbTimeSteps = len(outputs)
    time=np.zeros(nbTimeSteps)
    x = np.zeros([nbTimeSteps, 3])
    for i in range(0,nbTimeSteps):
        time[i] = outputs[i][0]
        x[i,:] = outputs[i][1]

    fig2 = plt.figure(figsize=(16, 10)) # in inches!
    plt.plot(time,x[:,0], 'red', label='Susceptible')
    plt.plot(time,x[:,1], 'green', label='Infected')
    plt.plot(time,x[:,2], 'blue', label='Recovered')
    plt.legend()
    plt.xlabel('Arbitrary time unit')
    plt.ylabel('Population numbers')
    string_for_title='SIRHigherOrderWellMixed, k_SI ='+str(params['k_SI']) +', k_IR = '+ str(params['k_IR'])\
      + ', S_ini = '+ str(params['S_ini']) + ', I_ini = '+str( params['I_ini']) + ', R_ini = '+ str(params['R_ini'])#\
      # + ', n_order = '+str(params['n_order'])
    plt.title(string_for_title)
    string_for_savefig = 'SIRWellMixed_k_SI='+str(params['k_SI']) +'_k_IR='+ str(params['k_IR'])\
      + '_S_ini='+ str(params['S_ini']) + '_I_ini='+str( params['I_ini']) + '_R_ini='+ str(params['R_ini'])\
      + 'n_order='+str(params['n_order'])
    plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
    # plt.close()
    plt.show()
def plot_SIRHigherOrderSeveralReactionsWellMixed(params,outputs):
    """
    plot SIRHigherOrderSeveralReactionsWellMixed model results
    S,I, R populations over time
    """
    # save outputs
    f = open(params['outputs_dir']+'/outputs.txt',"w")
    f.write(str(outputs))
    f.close()
    
    # save params inputs
    f = open(params['outputs_dir']+'/params.txt',"w")
    f.write(str(params))
    f.close()
    
    nbTimeSteps = len(outputs)
    time=np.zeros(nbTimeSteps)
    x = np.zeros([nbTimeSteps, 3])
    for i in range(0,nbTimeSteps):
        time[i] = outputs[i][0]
        x[i,:] = outputs[i][1]

    fig2 = plt.figure(figsize=(16, 10)) # in inches!
    plt.plot(time,x[:,0], 'red', label='Susceptible')
    plt.plot(time,x[:,1], 'green', label='Infected')
    plt.plot(time,x[:,2], 'blue', label='Recovered')
    plt.legend()
    plt.xlabel('Arbitrary time unit')
    plt.ylabel('Population numbers')
    string_for_title='SIRHigherOrderSeveralReactionsWellMixed, tilde_k_SI ='+str(params['tilde_k_SI']) + ', n_order = '+str(params['n_order'])      
    for i in range(0,params["n_order"]+1):
        str_to_search = 'p_'+str(i)
        p = params[str_to_search]
        string_for_title += ', '+str(str_to_search)+'='+str(p)
    plt.title(string_for_title)
    string_for_savefig = 'SIRHigherOrderSeveralReactionsWellMixed_tilde_k_SI='+str(params['tilde_k_SI']) + '_n_order='+str(params['n_order'])
    for i in range(0,params["n_order"]+1):
        str_to_search = 'p_'+str(i)
        p = params[str_to_search]
        string_for_savefig += '_'+str(str_to_search)+'='+str(p)
    plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
    # plt.close()
    plt.show()
def plot_SIHigherOrderSeveralReactionsWellMixed(params,outputs):
    """
    plot SIHigherOrderSeveralReactionsWellMixed model results
    S,I, R populations over time
    """
    # save outputs
    f = open(params['outputs_dir']+'/outputs.txt',"w")
    f.write(str(outputs))
    f.close()
    
    # save params inputs
    f = open(params['outputs_dir']+'/params.txt',"w")
    f.write(str(params))
    f.close()
    
    nbTimeSteps = len(outputs)
    time=np.zeros(nbTimeSteps)
    x = np.zeros([nbTimeSteps, 2])
    for i in range(0,nbTimeSteps):
        time[i] = outputs[i][0]
        x[i,:] = outputs[i][1]

    fig2 = plt.figure(figsize=(16, 10)) # in inches!
    plt.plot(time,x[:,0], 'red', label='Susceptible')
    plt.plot(time,x[:,1], 'green', label='Infected')    
    plt.legend()
    plt.xlabel('Arbitrary time unit')
    plt.ylabel('Population numbers')
    string_for_title='SIHigherOrderSeveralReactionsWellMixed, tilde_k_SI ='+str(params['tilde_k_SI']) + ', n_order = '+str(params['n_order'])
    for i in range(0,params["n_order"]+1):
        str_to_search = 'p_'+str(i)
        p = params[str_to_search]
        string_for_title += ', '+str(str_to_search)+'='+str(p)
    plt.title(string_for_title)
    string_for_savefig = 'SIHigherOrderSeveralReactionsWellMixed_tilde_k_SI='+str(params['tilde_k_SI']) + '_n_order='+str(params['n_order'])
    for i in range(0,params["n_order"]+1):
        str_to_search = 'p_'+str(i)
        p = params[str_to_search]
        string_for_savefig += '_'+str(str_to_search)+'='+str(p)
    plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
    # plt.close()
    plt.show()
def plot_SIRDWellMixed(params,outputs):
    """
    plot SIRWellMixed model results
    S,I, R, D populations over time
    """
    # save outputs
    f = open(params['outputs_dir']+'/outputs.txt',"w")
    f.write(str(outputs))
    f.close()
    
    # save params inputs
    f = open(params['outputs_dir']+'/params.txt',"w")
    f.write(str(params))
    f.close()
    
    nbTimeSteps = len(outputs)
    time=np.zeros(nbTimeSteps)
    x = np.zeros([nbTimeSteps, 4])
    
    for i in range(0,nbTimeSteps):
        time[i] = outputs[i][0]
        x[i,:] = outputs[i][1]
        
    fig2 = plt.figure(figsize=(16, 10)) # in inches!
    plt.plot(time,x[:,0], 'red', label='Susceptible')
    plt.plot(time,x[:,1], 'green', label='Infected')
    plt.plot(time,x[:,2], 'blue', label='Recovered')
    plt.plot(time,x[:,3], 'black', label='Dead')
    plt.legend()
    plt.xlabel('Arbitrary time unit')
    plt.ylabel('Population numbers')
    string_for_title='SIRDWellMixed, k_SI ='+str(params['k_SI']) +', k_IR = '+ str(params['k_IR'])+', k_ID = '+ str(params['k_ID']) + ', S_ini = '+ str(params['S_ini']) + ', I_ini = '+str( params['I_ini']) + ', R_ini = '+ str(params['R_ini'])+ ', D_ini = '+ str(params['D_ini'])
    plt.title(string_for_title)
    string_for_savefig = 'SIRDWellMixed_k_SI='+str(params['k_SI']) +'_k_IR='+ str(params['k_IR'])+'_k_ID='+ str(params['k_ID']) + '_S_ini='+ str(params['S_ini']) + '_I_ini='+str( params['I_ini']) + '_R_ini='+ str(params['R_ini'])+ '_D_ini='+ str(params['D_ini'])
    plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')    
def plot_SIRerlangWellMixed(params,outputs):
    """
    plot SIRerlangWellMixed model results
    S,I, R populations over time
    respective indices for plot are 0, m, nbStates
    """
    # save outputs
    f = open(params['outputs_dir']+'/outputs.txt',"w")
    f.write(str(outputs))
    f.close()
    
    # save params inputs
    f = open(params['outputs_dir']+'/params.txt',"w")
    f.write(str(params))
    f.close()
    
    nbTimeSteps = len(outputs)
    time=np.zeros(nbTimeSteps)
    nbStates = outputs[0][1].shape[0]
    m = params["nb_intermediate_states"]
    x = np.zeros([nbTimeSteps, nbStates])
    fig2 = plt.figure(figsize=(16, 10)) # in inches!
    for i in range(0,nbTimeSteps):
        time[i] = outputs[i][0]
        x[i,:] = outputs[i][1]        
    
    for ii in range(0,params['nbStates']):
        if (ii ==0):
            _label='Susceptible'
        if (ii==m):
            _label='Infected'
        if(ii==nbStates  - 1):
            _label='Recovered'
        if((ii>0) and (ii<m)):
            _label='I'+str(ii)
        if ((ii>m) and (ii<nbStates -1)):
            _label='R'+str(ii-m)        
        plt.plot(time,x[:,ii], label=_label)        
        plt.legend()
        plt.xlabel('Arbitrary time unit')
        plt.ylabel('Population numbers')
    string_for_title='SIRerlangWellMixed, k_SI ='+str(params['k_SI']) +', k_IR = '+ str(params['k_IR'])+' , nb_intermediate_states = ' +str(params['nb_intermediate_states']) + ', S_ini = '+ str(params['S_ini']) + ', I_ini = '+str( params['I_ini']) + ', R_ini = '+ str(params['R_ini'])
    plt.title(string_for_title)
    string_for_savefig = 'SIRerlangWellMixed_k_SI='+str(params['k_SI']) +'_k_IR='+ str(params['k_IR'])+'_nb_intermediate_states=' + str(params['nb_intermediate_states']) + '_S_ini='+ str(params['S_ini']) + '_I_ini='+str( params['I_ini']) + '_R_ini='+ str(params['R_ini'])
    plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
    plt.show()
def plot_SIRerlang1WellMixed(params,outputs):
    """
    plot SIRerlangWellMixed model results
    S,I, R populations over time
    respective indices for plot are 0, m, nbStates
    """
    # save outputs
    f = open(params['outputs_dir']+'/outputs.txt',"w")
    f.write(str(outputs))
    f.close()
    
    # save params inputs
    f = open(params['outputs_dir']+'/params.txt',"w")
    f.write(str(params))
    f.close()
    
    nbTimeSteps = len(outputs)
    time=np.zeros(nbTimeSteps)
    nbStates = outputs[0][1].shape[0]
    m = params["nb_intermediate_states"]
    x = np.zeros([nbTimeSteps, nbStates])
    fig2 = plt.figure(figsize=(16, 10)) # in inches!
    for i in range(0,nbTimeSteps):
        time[i] = outputs[i][0]
        x[i,:] = outputs[i][1]        
    
    tmp_sum = np.zeros(nbTimeSteps)
    for ii in range(0,params['nbStates']):
        if (ii ==0):
            tmp_sum+=x[:,ii]
            _label='Susceptible'
        if (ii==m):
            _label='Infected'
            tmp_sum+=x[:,ii]
        if(ii==nbStates  - 1):
            _label='Recovered'            
        if((ii>0) and (ii<m)):
            _label='I'+str(ii)     
            tmp_sum+=x[:,ii]
        plt.plot(time,x[:,ii], label=_label)            
    
    plt.plot(time,tmp_sum,color='black',label='sum of S + I_i (i=1:m-1)')    
    plt.legend()
    plt.xlabel('Arbitrary time unit')
    plt.ylabel('Population numbers')
    string_for_title='SIRerlang1WellMixed, k_SI ='+str(params['k_SI']) +', k_IR = '+ str(params['k_IR'])+' , nb_intermediate_states = ' +str(params['nb_intermediate_states']) + ', S_ini = '+ str(params['S_ini']) + ', I_ini = '+str( params['I_ini']) + ', R_ini = '+ str(params['R_ini'])
    plt.title(string_for_title)
    string_for_savefig = 'SIRerlang1WellMixed_k_SI='+str(params['k_SI']) +'_k_IR='+ str(params['k_IR'])+'_nb_intermediate_states=' + str(params['nb_intermediate_states']) + '_S_ini='+ str(params['S_ini']) + '_I_ini='+str( params['I_ini']) + '_R_ini='+ str(params['R_ini'])
    plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
    plt.show()  

def plot_SIRerlang2WellMixed(params,outputs):
    """
    plot SIRerlangWellMixed model results
    S,I, R populations over time
    respective indices for plot are 0, m, nbStates
    """
    # save outputs
    f = open(params['outputs_dir']+'/outputs.txt',"w")
    f.write(str(outputs))
    f.close()
    
    # save params inputs
    f = open(params['outputs_dir']+'/params.txt',"w")
    f.write(str(params))
    f.close()
    
    nbTimeSteps = len(outputs)
    time=np.zeros(nbTimeSteps)
    nbStates = outputs[0][1].shape[0]
    m = params["nb_intermediate_states"]
    x = np.zeros([nbTimeSteps, nbStates])
    fig2 = plt.figure(figsize=(16, 10)) # in inches!
    for i in range(0,nbTimeSteps):
        time[i] = outputs[i][0]
        x[i,:] = outputs[i][1]        
    
    for ii in range(0,params['nbStates']):
        if (ii ==0):
            _label='Susceptible'
        if (ii==1):
            _label='Infected'
        if(ii==nbStates  - 1):
            _label='Recovered'
        if((ii>1) and (ii<params['nbStates'] - 1)):
            _label='R'+str(ii-1)     
        plt.plot(time,x[:,ii], label=_label)        
        plt.legend()
        plt.xlabel('Arbitrary time unit')
        plt.ylabel('Population numbers')
    string_for_title='SIRerlang2WellMixed, k_SI ='+str(params['k_SI']) +', k_IR = '+ str(params['k_IR'])+' , nb_intermediate_states = ' +str(params['nb_intermediate_states']) + ', S_ini = '+ str(params['S_ini']) + ', I_ini = '+str( params['I_ini']) + ', R_ini = '+ str(params['R_ini'])
    plt.title(string_for_title)
    string_for_savefig = 'SIRerlang2WellMixed_k_SI='+str(params['k_SI']) +'_k_IR='+ str(params['k_IR'])+'_nb_intermediate_states=' + str(params['nb_intermediate_states']) + '_S_ini='+ str(params['S_ini']) + '_I_ini='+str( params['I_ini']) + '_R_ini='+ str(params['R_ini'])
    plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
    plt.close()
    # plt.show()
    
    if ('WaitingTimes' in params['plot_options']):
        count_I_m = 1
        count_S = params['S_ini']
        count_R = 0
        S = x[:,0]
        I_m = x[:,nbStates-2]
        R = x[:,nbStates-1]
        arrival_times = []
        count_append = 0
        for i in range(0,nbTimeSteps):
            if (R[i] > count_R):
                count_R +=1
                arrival_times.append(time[i])
        
        # convert to numpy
        arrival_times = np.array(arrival_times)
        
        # m * rate gives rate of erlang distrib
        k_IR = params['k_IR'] 
        m = params['nb_intermediate_states'] 
        Lambda_IR = m * k_IR
        time_real = np.linspace(0, np.max(arrival_times), 1000)        
        # compute (m-1)!
        facto_m_minus_1 = 1
        for i in range(1,m):
            facto_m_minus_1 *= i
        real_erlang_distrib = Lambda_IR * (Lambda_IR * time_real) **(m-1) * np.exp(-Lambda_IR * time_real) / facto_m_minus_1        
        
        # mu = 38.3479
        # sigma = 6.5
        # gaussian_distrib = 1./(sigma * np.sqrt(2 * np.pi)) * np.exp(-0.5*(time_real-mu)**2 / sigma**2)
        
        # plot histogram
        fig100 = plt.figure(figsize=(16, 10)) # in inches!
        plt.hist(arrival_times, bins=20, density=True, label='Waiting times distribution')
        # plt.hist(departure_times, bins=int(params['S_ini']/200), density=True, label='Waiting times distribution (departure)')
        plt.xlabel('Waiting times')
        plt.ylabel('PDF')
        # plot real erlang distrib
        plt.plot(time_real, real_erlang_distrib, label='Erlang distribution')
        # plt.plot(time_real, gaussian_distrib, label='Gaussian distribution')
        plt.legend()
        string_for_title='SIRerlang2WellMixed k_SI ='+str(params['k_SI'])+' k_IR = '+str(params['k_IR'])+' , nb_intermediate_states = ' +str(params['nb_intermediate_states']) + ', S_ini = '+ str(params['S_ini']) + ', I_ini = '+str( params['I_ini'])
        plt.title(string_for_title)
        string_for_savefig = 'SIRerlang2WellMixed_waiting_times_pdf_k_SI='+str(params['k_SI']) +'_k_IR='+str(params['k_IR'])+'_nb_intermediate_states=' + str(params['nb_intermediate_states']) + '_S_ini='+ str(params['S_ini']) + '_I_ini='+str( params['I_ini'])
        plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
        # plt.close()
        plt.show()   
      
def plot_SIWellMixed(params,outputs):
    """
    plot SIRerlangWellMixed model results
    S,I
    respective indices for plot are 0, m, nbStates
    """
    # save outputs
    f = open(params['outputs_dir']+'/outputs.txt',"w")
    f.write(str(outputs))
    f.close()
    
    # save params inputs
    f = open(params['outputs_dir']+'/params.txt',"w")
    f.write(str(params))
    f.close()
    
    nbTimeSteps = len(outputs)
    time=np.zeros(nbTimeSteps)
    nbStates = outputs[0][1].shape[0]
    x = np.zeros([nbTimeSteps, nbStates])
    fig2 = plt.figure(figsize=(16, 10)) # in inches!
    for i in range(0,nbTimeSteps):
        time[i] = outputs[i][0]
        x[i,:] = outputs[i][1]        
    
    I_sum = np.zeros(nbTimeSteps)
    for ii in range(0,params['nbStates']):
        if (ii ==0):
            _label='Susceptible'
        if (ii==1):
            _label='Infected'                    
        plt.plot(time,x[:,ii], label=_label)        
    plt.legend()
    plt.xlabel('Arbitrary time unit')
    plt.ylabel('Population numbers')
    string_for_title='SIWellMixed, k_SI ='+str(params['k_SI'])+ ', S_ini = '+ str(params['S_ini']) + ', I_ini = '+str( params['I_ini']) 
    plt.title(string_for_title)
    string_for_savefig = 'SIWellMixed_k_SI='+str(params['k_SI'])+ '_S_ini='+ str(params['S_ini']) + '_I_ini='+str( params['I_ini'])
    plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
    plt.show()   

def plot_SIAndIntermediateErlangWellMixed(params,outputs):
    """
    plot SIAndIntermediateErlangWellMixed model results
    S,I
    respective indices for plot are 0, m, nbStates
    """
    # save outputs
    fname_outputs  = params['outputs_dir']+'/outputs_'+params['model_name']+'_m='+str(params['nb_intermediate_states'])+'.txt'
    f = open(fname_outputs,"w")
    f.write(str(outputs))
    f.close()
    
    # save params inputs
    fname_params  = params['outputs_dir']+'/params_'+params['model_name']+'_m='+str(params['nb_intermediate_states'])+'.txt'
    f = open(fname_params,"w")
    f.write(str(params))
    f.close()
    
    nbTimeSteps = len(outputs)
    time=np.zeros(nbTimeSteps)
    nbStates = outputs[0][1].shape[0]
    m = params["nb_intermediate_states"]
    x = np.zeros([nbTimeSteps, nbStates])
    fig2 = plt.figure(figsize=(16, 10)) # in inches!
    for i in range(0,nbTimeSteps):
        time[i] = outputs[i][0]
        x[i,:] = outputs[i][1]
    
    tmp_sum = np.zeros(nbTimeSteps)
    for ii in range(0,params['nbStates']):
        if (ii ==0):
            _label='Susceptible'
            tmp_sum += x[:,ii]
        if (ii ==1):
            _label='Intermediate L'
            tmp_sum += x[:,ii]
        if (ii==params['nbStates']-1):
            _label='Infected'     
        if((ii>1) and (ii<params['nbStates']-1)):
            _label='I'+str(ii-1)     
            tmp_sum += x[:,ii]
        plt.plot(time,x[:,ii], label=_label)            
    
    plt.plot(time,tmp_sum,color='black',label='sum of S + L + I_i (i=1:m-1)')
    plt.legend()
    plt.xlabel('Arbitrary time unit')
    plt.ylabel('Population numbers')
    string_for_title='SIAndIntermediateErlangWellMixed, k_SL ='+str(params['k_SL'])+' k_LI ='+str(params['k_LI'])+' , nb_intermediate_states = ' +str(params['nb_intermediate_states']) + ', S_ini = '+ str(params['S_ini']) + ', I_ini = '+str( params['I_ini'])
    plt.title(string_for_title)
    string_for_savefig = 'SIAndIntermediateErlangWellMixed_k_SL='+str(params['k_SL'])+str(params['k_LI'])+'_k_LI=' +'_nb_intermediate_states=' + str(params['nb_intermediate_states']) + '_S_ini='+ str(params['S_ini']) + '_I_ini='+str( params['I_ini'])
    plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
    # plt.close()
    plt.show()  
    
    # plot waiting times    
    if ('WaitingTimes' in params['plot_options']):
        # FOR NOW: it only works in the case where S = 0, I_m = 0, L = L_ini
        # if (params['S_ini'] == 0) and (params['I_ini'] == 0):
        count_I_m = 0
        I_m = x[:,nbStates-1]
        arrival_times = []
        for i in range(0,nbTimeSteps):
            if (I_m[i] > count_I_m) :
                count_I_m +=1
                arrival_times.append(time[i])
        # convert to numpy
        arrival_times = np.array(arrival_times)        
        
        # m * rate gives rate of erlang distrib
        k_LI = params['k_LI']
        k_SL = params['k_SL']
        S_ini = params['S_ini']
        m = params['nb_intermediate_states']
        # Lambda_LI = m * k_LI
        Lambda_LI = m * k_SL
        Lambda_LI = m * k_SL * S_ini
        time_real = np.linspace(0, np.max(arrival_times), 1000)        
        # compute (m-1)!
        facto_m_minus_1 = 1
        for i in range(1,m):
            facto_m_minus_1 *= i
        real_erlang_distrib = Lambda_LI * (Lambda_LI * time_real) **(m-1) * np.exp(-Lambda_LI * time_real) / facto_m_minus_1
        
        # plot histogram
        fig100 = plt.figure(figsize=(16, 10)) # in inches!
        plt.hist(arrival_times, bins=int(200), density=True, label='Waiting times distribution')
        plt.xlabel('Waiting times')
        plt.ylabel('PDF')
        # plot real erlang distrib
        plt.plot(time_real, real_erlang_distrib, label='Erlang distribution')
        plt.legend()
        string_for_title='SIAndIntermediateErlangWellMixed, k_SL ='+str(params['k_SL'])+' k_LI ='+str(params['k_LI'])+' , nb_intermediate_states = ' +str(params['nb_intermediate_states']) + ', S_ini = '+ str(params['S_ini']) + ', I_ini = '+str( params['I_ini']) + ', L_ini = '+str( params['L_ini'])
        plt.title(string_for_title)
        string_for_savefig = 'SIAndIntermediateErlangWellMixed_waiting_times_pdf_k_SL='+str(params['k_SL'])+str(params['k_LI'])+'_k_LI=' +'_nb_intermediate_states=' + str(params['nb_intermediate_states']) + '_S_ini='+ str(params['S_ini']) + '_I_ini='+str( params['I_ini'])+ '_L_ini='+str( params['L_ini'])
        plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
        # plt.close()
        plt.show()
              
def plot_SIerlangWellMixed(params,outputs):
    """
    plot SIerlangWellMixed model results
    S,I
    respective indices for plot are 0, m, nbStates
    """
    # save outputs
    fname_outputs  = params['outputs_dir']+'/outputs_'+params['model_name']+'_m='+str(params['nb_intermediate_states'])+'.txt'
    f = open(fname_outputs,"w")
    f.write(str(outputs))
    f.close()
    
    # save params inputs
    fname_params  = params['outputs_dir']+'/params_'+params['model_name']+'_m='+str(params['nb_intermediate_states'])+'.txt'
    f = open(fname_params,"w")
    f.write(str(params))
    f.close()
    
    nbTimeSteps = len(outputs)
    time=np.zeros(nbTimeSteps)
    nbStates = outputs[0][1].shape[0]
    m = params["nb_intermediate_states"]
    x = np.zeros([nbTimeSteps, nbStates])
    fig2 = plt.figure(figsize=(16, 10)) # in inches!
    for i in range(0,nbTimeSteps):
        time[i] = outputs[i][0]
        x[i,:] = outputs[i][1]        
    
    tmp_sum = np.zeros(nbTimeSteps)
    for ii in range(0,params['nbStates']):
        if (ii ==0):
            _label='Susceptible'
            tmp_sum += x[:,ii]
        if (ii==m):
            _label='Infected'            
        if((ii>0) and (ii<m)):
            _label='I'+str(ii)     
            tmp_sum += x[:,ii]
        plt.plot(time,x[:,ii], label=_label)            
    
    plt.plot(time,tmp_sum,color='black',label='sum of S + I_i (i=1:m-1)')
    plt.legend()
    plt.xlabel('Arbitrary time unit')
    plt.ylabel('Population numbers')
    # plt.yscale('log')
    string_for_title='SIRerlangWellMixed, k_SI ='+str(params['k_SI'])+' , nb_intermediate_states = ' +str(params['nb_intermediate_states']) + ', S_ini = '+ str(params['S_ini']) + ', I_ini = '+str( params['I_ini'])
    plt.title(string_for_title)
    string_for_savefig = 'SIRerlangWellMixed_k_SI='+str(params['k_SI']) +'_nb_intermediate_states=' + str(params['nb_intermediate_states']) + '_S_ini='+ str(params['S_ini']) + '_I_ini='+str( params['I_ini'])
    plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
    # plt.close() 
    plt.show()  
    # plot waiting times    
    if ('WaitingTimes' in params['plot_options']):
        count_I_m = 1
        I_m = x[:,nbStates-1]
        arrival_times = []
        for i in range(0,nbTimeSteps):
            if (I_m[i] > count_I_m):
                arrival_times.append(time[i])
                count_I_m +=1
                
        # convert to numpy        
        arrival_times =  np.array(arrival_times)
            
        # m * rate gives rate of erlang distrib
        k_SI = params['k_SI'] 
        m = params['nb_intermediate_states']
        S_ini = params['S_ini']
        # Lambda_SI = m * k_SI * S_ini
        Lambda_SI = m * k_SI * S_ini
        time_real = np.linspace(0, np.max(arrival_times), 1000)        
        # compute (m-1)!
        facto_m_minus_1 = 1
        for i in range(1,m):
            facto_m_minus_1 *= i
        real_erlang_distrib = Lambda_SI * (Lambda_SI * time_real) **(m-1) * np.exp(-Lambda_SI * time_real) / facto_m_minus_1        
        
        # plot histogram
        fig100 = plt.figure(figsize=(16, 10)) # in inches!
        plt.hist(arrival_times, bins=int(200), density=True, label='Waiting times distribution (arrival)')
        plt.xlabel('Waiting times')
        plt.ylabel('PDF')
        # plot real erlang distrib
        plt.plot(time_real, real_erlang_distrib, label='Erlang distribution')
        plt.legend()
        string_for_title='SIerlangWellMixed k_SI ='+str(params['k_SI'])+' , nb_intermediate_states = ' +str(params['nb_intermediate_states']) + ', S_ini = '+ str(params['S_ini']) + ', I_ini = '+str( params['I_ini'])
        plt.title(string_for_title)
        string_for_savefig = 'SIerlangWellMixed_waiting_times_pdf_k_SI='+str(params['k_SI']) +'_nb_intermediate_states=' + str(params['nb_intermediate_states']) + '_S_ini='+ str(params['S_ini']) + '_I_ini='+str( params['I_ini'])
        plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
        # plt.close()
        plt.show()   

def plot_SIerlangWellMixedModified(params,outputs, waiting_times):
    """
    plot SIerlangWellMixed model results
    S,I
    respective indices for plot are 0, m, nbStates
    """
    # save outputs
    fname_outputs  = params['outputs_dir']+'/outputs_'+params['model_name']+'_m='+str(params['nb_intermediate_states'])+'.txt'
    f = open(fname_outputs,"w")
    f.write(str(outputs))
    f.close()
    
    # save params inputs
    fname_params  = params['outputs_dir']+'/params_'+params['model_name']+'_m='+str(params['nb_intermediate_states'])+'.txt'
    f = open(fname_params,"w")
    f.write(str(params))
    f.close()
    
    nbTimeSteps = len(outputs)
    time=np.zeros(nbTimeSteps)
    nbStates = outputs[0][1].shape[0]
    m = params["nb_intermediate_states"]
    x = np.zeros([nbTimeSteps, nbStates])
    fig2 = plt.figure(figsize=(16, 10)) # in inches!
    for i in range(0,nbTimeSteps):
        time[i] = outputs[i][0]
        x[i,:] = outputs[i][1]        
    
    tmp_sum = np.zeros(nbTimeSteps)
    for ii in range(0,params['nbStates']):
        if (ii ==0):
            _label='Susceptible'
            tmp_sum += x[:,ii]
        if (ii==m):
            _label='Infected'            
        if((ii>0) and (ii<m)):
            _label='I'+str(ii)     
            tmp_sum += x[:,ii]
        plt.plot(time,x[:,ii], label=_label)            
    
    plt.plot(time,tmp_sum,color='black',label='sum of S + I_i (i=1:m-1)')
    plt.legend()
    plt.xlabel('Arbitrary time unit')
    plt.ylabel('Population numbers')
    # plt.yscale('log')
    string_for_title='SIRerlangWellMixedModified, k_SI ='+str(params['k_SI'])+' , nb_intermediate_states = ' +str(params['nb_intermediate_states']) + ', S_ini = '+ str(params['S_ini']) + ', I_ini = '+str( params['I_ini'])
    plt.title(string_for_title)
    string_for_savefig = 'SIRerlangWellMixedModified_k_SI='+str(params['k_SI']) +'_nb_intermediate_states=' + str(params['nb_intermediate_states']) + '_S_ini='+ str(params['S_ini']) + '_I_ini='+str( params['I_ini'])
    plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
    # plt.close() 
    plt.show()  
    # plot waiting times    
    if ('WaitingTimes' in params['plot_options']):
        # m * rate gives rate of erlang distrib
        k_SI = params['k_SI'] 
        m = params['nb_intermediate_states']
        S_ini = params['S_ini']
        I_ini = params['I_ini']
        Lambda_SI = (m) * k_SI
        time_real = np.linspace(0, np.max(waiting_times), 1000)
        # compute (m-1)!
        facto_m_minus_1 = 1
        for i in range(1,m):
            facto_m_minus_1 *= i
        real_erlang_distrib = Lambda_SI * (Lambda_SI * time_real) **(m-1) * np.exp(-Lambda_SI * time_real) / facto_m_minus_1        
                
        # plot histogram
        fig100 = plt.figure(figsize=(16, 10)) # in inches!
        # plt.hist(arrival_times, bins=int(params['S_ini']/200), density=True, label='Waiting times distribution (arrival)')
        # plt.hist(departure_times, bins=int(params['S_ini']/200), density=True, label='Waiting times distribution (departure)')
        plt.hist(waiting_times, bins=int(params['S_ini']/100), density=True, label='Waiting times distribution')
        plt.xlabel('Waiting times')
        plt.ylabel('PDF')
        # plot real erlang distrib
        plt.plot(time_real, real_erlang_distrib, label='Erlang distribution')
        # plt.plot(time_real, gaussian_distrib, label='Gaussian distribution')
        plt.legend()
        string_for_title='SIRerlangWellMixedModified waiting times k_SI ='+str(params['k_SI'])+' , nb_intermediate_states = ' +str(params['nb_intermediate_states']) + ', S_ini = '+ str(params['S_ini']) + ', I_ini = '+str( params['I_ini'])
        plt.title(string_for_title)
        string_for_savefig = 'SIRerlangWellMixedModified_waiting_times_pdf_k_SI='+str(params['k_SI']) +'_nb_intermediate_states=' + str(params['nb_intermediate_states']) + '_S_ini='+ str(params['S_ini']) + '_I_ini='+str( params['I_ini'])
        plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
        # plt.close()
        plt.show()   

def plot_SIAndIntermediateErlangWellMixedModified(params,outputs, waiting_times):
    """
    plot SIerlangWellMixed model results
    S,I
    respective indices for plot are 0, m, nbStates
    """
    # save outputs
    fname_outputs  = params['outputs_dir']+'/outputs_'+params['model_name']+'_m='+str(params['nb_intermediate_states'])+'.txt'
    f = open(fname_outputs,"w")
    f.write(str(outputs))
    f.close()
    
    # save params inputs
    fname_params  = params['outputs_dir']+'/params_'+params['model_name']+'_m='+str(params['nb_intermediate_states'])+'.txt'
    f = open(fname_params,"w")
    f.write(str(params))
    f.close()
    
    nbTimeSteps = len(outputs)
    time=np.zeros(nbTimeSteps)
    nbStates = outputs[0][1].shape[0]
    m = params["nb_intermediate_states"]
    x = np.zeros([nbTimeSteps, nbStates])
    fig2 = plt.figure(figsize=(16, 10)) # in inches!
    for i in range(0,nbTimeSteps):
        time[i] = outputs[i][0]
        x[i,:] = outputs[i][1]        
    
    tmp_sum = np.zeros(nbTimeSteps)
    for ii in range(0,params['nbStates']):
        if (ii ==0):
            _label='Susceptible'
            tmp_sum += x[:,ii]
        if (ii==m):
            _label='Infected'            
        if((ii>0) and (ii<m)):
            _label='I'+str(ii)     
            tmp_sum += x[:,ii]
        plt.plot(time,x[:,ii], label=_label)            
    
    plt.plot(time,tmp_sum,color='black',label='sum of S + I_i (i=1:m-1)')
    plt.legend()
    plt.xlabel('Arbitrary time unit')
    plt.ylabel('Population numbers')
    # plt.yscale('log')
    string_for_title='SIAndIntermediateErlangWellMixedModified '+str(params['nb_intermediate_states']) + ', S_ini = '+ str(params['S_ini'])+ ', L_ini = '+str( params['L_ini']) + ', I_ini = '+str( params['I_ini'])
    plt.title(string_for_title)
    string_for_savefig = 'SIAndIntermediateErlangWellMixedModified_k_SL='+str(params['k_SL'])+'_k_LI='+str(params['k_LI']) +'_nb_intermediate_states=' + str(params['nb_intermediate_states']) + '_S_ini='+ str(params['S_ini'])+ '_L_ini='+ str(params['L_ini']) + '_I_ini='+str( params['I_ini'])
    plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
    # plt.close() 
    plt.show()  
    # plot waiting times    
    if ('WaitingTimes' in params['plot_options']):
        # m * rate gives rate of erlang distrib
        k_SL = params['k_SL'] 
        k_LI = params['k_LI'] 
        m = params['nb_intermediate_states']
        S_ini = params['S_ini']
        I_ini = params['I_ini']
        Lambda_LI = m*k_SL
        time_real = np.linspace(0, np.max(waiting_times), 1000)
        # compute (m-1)!
        facto_m_minus_1 = 1
        for i in range(1,m):
            facto_m_minus_1 *= i
        real_erlang_distrib = Lambda_LI * (Lambda_LI * time_real) **(m-1) * np.exp(-Lambda_LI * time_real) / facto_m_minus_1        
                
        # plot histogram
        fig100 = plt.figure(figsize=(16, 10)) # in inches!
        plt.hist(waiting_times, bins=int(params['S_ini']/100), density=True, label='Waiting times distribution')
        plt.xlabel('Waiting times')
        plt.ylabel('PDF')
        # plot real erlang distrib
        plt.plot(time_real, real_erlang_distrib, label='Erlang distribution')
        plt.legend()
        string_for_title='SIRerlangWellMixedModified waiting times, nb_intermediate_states = ' +str(params['nb_intermediate_states']) + ', S_ini = '+ str(params['S_ini']) + ', I_ini = '+str( params['I_ini'])
        plt.title(string_for_title)
        string_for_savefig = 'SIRerlangWellMixedModified_waiting_times_pdf_k_SL='+str(params['k_SL'])+'_k_LI='+str(params['k_LI']) +'_nb_intermediate_states=' + str(params['nb_intermediate_states']) + '_S_ini='+ str(params['S_ini'])+ '_L_ini='+str( params['L_ini']) + '_I_ini='+str( params['I_ini'])
        plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
        # plt.close()
        plt.show()   

def plot_SIerlangFixSWellMixed(params,outputs):
    """
    plot SIerlangWellMixed model results
    S,I
    respective indices for plot are 0, m, nbStates
    """
    # save outputs
    fname_outputs  = params['outputs_dir']+'/outputs_'+params['model_name']+'_m='+str(params['nb_intermediate_states'])+'.txt'
    f = open(fname_outputs,"w")
    f.write(str(outputs))
    f.close()
    
    # save params inputs
    fname_params  = params['outputs_dir']+'/params_'+params['model_name']+'_m='+str(params['nb_intermediate_states'])+'.txt'
    f = open(fname_params,"w")
    f.write(str(params))
    f.close()
    
    nbTimeSteps = len(outputs)
    time=np.zeros(nbTimeSteps)
    nbStates = outputs[0][1].shape[0]
    m = params["nb_intermediate_states"]
    x = np.zeros([nbTimeSteps, nbStates])
    fig2 = plt.figure(figsize=(16, 10)) # in inches!
    for i in range(0,nbTimeSteps):
        time[i] = outputs[i][0]
        x[i,:] = outputs[i][1]        
    
    tmp_sum = np.zeros(nbTimeSteps)
    for ii in range(0,params['nbStates']):
        if (ii ==0):
            _label='Susceptible'
            tmp_sum += x[:,ii]
        if (ii==m):
            _label='Infected'            
        if((ii>0) and (ii<m)):
            _label='I'+str(ii)     
            tmp_sum += x[:,ii]
        plt.plot(time,x[:,ii], label=_label)            
    
    plt.plot(time,tmp_sum,color='black',label='sum of S + I_i (i=1:m-1)')
    plt.legend()
    plt.xlabel('Arbitrary time unit')
    plt.ylabel('Population numbers')
    # plt.yscale('log')
    string_for_title='SIRerlangWellMixed, k_SI ='+str(params['k_SI'])+' , nb_intermediate_states = ' +str(params['nb_intermediate_states']) + ', S_ini = '+ str(params['S_ini']) + ', I_ini = '+str( params['I_ini'])
    plt.title(string_for_title)
    string_for_savefig = 'SIRerlangWellMixed_k_SI='+str(params['k_SI']) +'_nb_intermediate_states=' + str(params['nb_intermediate_states']) + '_S_ini='+ str(params['S_ini']) + '_I_ini='+str( params['I_ini'])
    plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
    # plt.close() 
    plt.show()  
    # plot waiting times    
    if ('WaitingTimes' in params['plot_options']):
        count_I_m = 1
        I_m = x[:,nbStates-1]
        arrival_times = []
        count_append = 0
        count_loop = 0
        for i in range(0,nbTimeSteps):
            # if (I_m[i] > count_I_m):
                # if (count_loop == 0):
                    # arrival_times.append(time[i])
                # else:
                    # arrival_times.append(time[i] - arrival_times[count_loop -1])

                arrival_times.append(time[i])
                count_I_m +=1
                count_loop += 1            
        # convert to numpy        
        arrival_times =  np.array(arrival_times)
            
        # m * rate gives rate of erlang distrib
        k_SI = params['k_SI'] 
        m = params['nb_intermediate_states'] 
        S_ini = params['S_ini']
        # Lambda_SI = m * k_SI * S_ini
        Lambda_SI = m * k_SI  * S_ini
        time_real = np.linspace(0, np.max(arrival_times), 1000)        
        # compute (m-1)!
        facto_m_minus_1 = 1
        for i in range(1,m):
            facto_m_minus_1 *= i
        real_erlang_distrib = Lambda_SI * (Lambda_SI * time_real) **(m-1) * np.exp(-Lambda_SI * time_real) / facto_m_minus_1        

        # plot histogram
        fig100 = plt.figure(figsize=(16, 10)) # in inches!
        plt.hist(arrival_times, bins=int(params['S_ini']/200), density=True, label='Waiting times distribution (arrival)')
        plt.xlabel('Waiting times')
        plt.ylabel('PDF')
        # plot real erlang distrib
        plt.plot(time_real, real_erlang_distrib, label='Erlang distribution')
        plt.legend()
        string_for_title='SIerlangFixSWellMixed k_SI ='+str(params['k_SI'])+' , nb_intermediate_states = ' +str(params['nb_intermediate_states']) + ', S_ini = '+ str(params['S_ini']) + ', I_ini = '+str( params['I_ini'])
        plt.title(string_for_title)
        string_for_savefig = 'SIerlangFixSWellMixed_waiting_times_pdf_k_SI='+str(params['k_SI']) +'_nb_intermediate_states=' + str(params['nb_intermediate_states']) + '_S_ini='+ str(params['S_ini']) + '_I_ini='+str( params['I_ini'])
        plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
        # plt.close()
        plt.show()   

def plot_SIerlangOtherWellMixed(params,outputs):
    """
    plot_SIerlangOtherWellMixed model results
    S,I
    respective indices for plot are 0, m, nbStates
    """
    # save outputs
    fname_outputs  = params['outputs_dir']+'/outputs_'+params['model_name']+'_m='+str(params['nb_intermediate_states'])+'.txt'
    f = open(fname_outputs,"w")
    f.write(str(outputs))
    f.close()
    
    # save params inputs
    fname_params  = params['outputs_dir']+'/params_'+params['model_name']+'_m='+str(params['nb_intermediate_states'])+'.txt'
    f = open(fname_params,"w")
    f.write(str(params))
    f.close()
    
    nbTimeSteps = len(outputs)
    time=np.zeros(nbTimeSteps)
    nbStates = outputs[0][1].shape[0]
    m = params["nb_intermediate_states"]
    x = np.zeros([nbTimeSteps, nbStates])
    fig2 = plt.figure(figsize=(16, 10)) # in inches!
    for i in range(0,nbTimeSteps):
        time[i] = outputs[i][0]
        x[i,:] = outputs[i][1]        
    
    tmp_sum = np.zeros(nbTimeSteps)
    for ii in range(0,params['nbStates']):
        if (ii ==0):
            _label='Susceptible'
            tmp_sum += x[:,ii]
        if (ii==m):
            _label='Infected'            
        if((ii>0) and (ii<m)):
            _label='I'+str(ii)     
            tmp_sum += x[:,ii]
        plt.plot(time,x[:,ii], label=_label)            
    
    plt.plot(time,tmp_sum,color='black',label='sum of S + I_i (i=1:m-1)')
    plt.legend()
    plt.xlabel('Arbitrary time unit')
    plt.ylabel('Population numbers')
    string_for_title='SIerlangOtherWellMixed, k_SI ='+str(params['k_SI'])+' , nb_intermediate_states = ' +str(params['nb_intermediate_states']) + ', S_ini = '+ str(params['S_ini']) + ', I_ini = '+str( params['I_ini'])
    plt.title(string_for_title)
    string_for_savefig = 'SIerlangOtherWellMixed_k_SI='+str(params['k_SI']) +'_nb_intermediate_states=' + str(params['nb_intermediate_states']) + '_S_ini='+ str(params['S_ini']) + '_I_ini='+str( params['I_ini'])
    plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
    plt.show()  
    
     # plot waiting times    
    if ('WaitingTimes' in params['plot_options']):
        count_I_m = 1
        I_m = x[:,nbStates-1]
        arrival_times = []
        for i in range(0,nbTimeSteps):
            if (I_m[i] > count_I_m):
                arrival_times.append(time[i])
                count_I_m +=1
                
        # convert to numpy        
        arrival_times =  np.array(arrival_times)
            
        # m * rate gives rate of erlang distrib
        k_SI = params['k_SI'] 
        m = params['nb_intermediate_states']
        S_ini = params['S_ini']
        # Lambda_SI = m * k_SI * S_ini
        Lambda_SI = m * k_SI #* S_ini
        time_real = np.linspace(0, np.max(arrival_times), 1000)        
        # compute (m-1)!
        facto_m_minus_1 = 1
        for i in range(1,m):
            facto_m_minus_1 *= i
        real_erlang_distrib = Lambda_SI * (Lambda_SI * time_real) **(m-1) * np.exp(-Lambda_SI * time_real) / facto_m_minus_1        
        
        # plot histogram
        fig100 = plt.figure(figsize=(16, 10)) # in inches!
        plt.hist(arrival_times, bins=int(40), density=True, label='Waiting times distribution (arrival)')
        plt.xlabel('Waiting times')
        plt.ylabel('PDF')
        # plot real erlang distrib
        plt.plot(time_real, real_erlang_distrib, label='Erlang distribution')
        plt.legend()
        string_for_title='SIerlangOtherWellMixed k_SI ='+str(params['k_SI'])+' , nb_intermediate_states = ' +str(params['nb_intermediate_states']) + ', S_ini = '+ str(params['S_ini']) + ', I_ini = '+str( params['I_ini'])
        plt.title(string_for_title)
        string_for_savefig = 'SIerlangOtherWellMixed_waiting_times_pdf_k_SI='+str(params['k_SI']) +'_nb_intermediate_states=' + str(params['nb_intermediate_states']) + '_S_ini='+ str(params['S_ini']) + '_I_ini='+str( params['I_ini'])
        plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
        # plt.close()
        plt.show()   

    
def plot_SILinearErlangWellMixed(params,outputs):
    """
    plot SIerlangWellMixed model results
    S,I
    respective indices for plot are 0, m, nbStates
    """
    # save outputs
    f = open(params['outputs_dir']+'/outputs.txt',"w")
    f.write(str(outputs))
    f.close()
    
    # save params inputs
    f = open(params['outputs_dir']+'/params.txt',"w")
    f.write(str(params))
    f.close()
    
    nbTimeSteps = len(outputs)
    time=np.zeros(nbTimeSteps)
    nbStates = outputs[0][1].shape[0]
    m = params["nb_intermediate_states"]
    x = np.zeros([nbTimeSteps, nbStates])
    fig2 = plt.figure(figsize=(16, 10)) # in inches!
    for i in range(0,nbTimeSteps):
        time[i] = outputs[i][0]
        x[i,:] = outputs[i][1]        
    
    tmp_sum = np.zeros(nbTimeSteps)
    for ii in range(0,params['nbStates']):
        if (ii ==0):
            _label='Susceptible'
            tmp_sum += x[:,ii]
        if (ii==m):
            _label='Infected'
        if((ii>0) and (ii<m)):
            _label='I'+str(ii)
            tmp_sum += x[:,ii]
        plt.plot(time,x[:,ii], label=_label)
    
    
    plt.plot(time,tmp_sum,color='gray',label='sum of S + I_i (i=1:m-1)')
    plt.plot(time,tmp_sum+x[:,m],color='black',label='sum of all')
    plt.legend()
    plt.xlabel('Arbitrary time unit')
    plt.ylabel('Population numbers')
    string_for_title='SIRerlang1WellMixed, k_SI ='+str(params['k_SI'])+' , nb_intermediate_states = ' +str(params['nb_intermediate_states']) + ', S_ini = '+ str(params['S_ini']) + ', I_ini = '+str( params['I_ini'])
    plt.title(string_for_title)
    string_for_savefig = 'SIRerlang1WellMixed_k_SI='+str(params['k_SI']) +'_nb_intermediate_states=' + str(params['nb_intermediate_states']) + '_S_ini='+ str(params['S_ini']) + '_I_ini='+str( params['I_ini'])
    plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
    plt.show()  

def plot_SIROnGraph(params,outputs, G):
    """
    plot plot_SIROnGraph model results
    S,I, R populations over time
    depending values contained in params['plot_options'] we can either:
    WellMixed: plot the well mixed limit by summing all compartments
    AllNodes: plot the dynamic on every node separately
    AllGraphs: plot all graphs over all time steps
    FinalGraph: plot only the final graph at the final time step
    """
    print('plot_SIROnGraph')
    # save outputs
    f = open(params['outputs_dir']+'/outputs.txt',"w")
    f.write(str(outputs))
    f.close()
    
    # save params inputs
    f = open(params['outputs_dir']+'/params.txt',"w")
    f.write(str(params))
    f.close()
    
    nbTimeSteps = len(outputs)
    time=np.zeros(nbTimeSteps)
    x = np.zeros([nbTimeSteps, 3])
    dic_graphs = {}
    for i in range(0,nbTimeSteps):
        time[i] = outputs[i][0]
        dic_x = outputs[i][1]        
        for n in dic_x.keys():
            x[i,0] += dic_x[n][0]
            x[i,1] += dic_x[n][1]
            x[i,2] += dic_x[n][2]            
    
    if ('WellMixed' in params['plot_options']):
        print("Running: plot option: 'WellMixed'")
        ### --- plot wellmixed model ( = sum over all compartments i.e. nodes of the graph)
        fig10 = plt.figure(figsize=(16, 10)) # in inches!
        plt.plot(time,x[:,0], 'red', label='Susceptible')
        plt.plot(time,x[:,1], 'green', label='Infected')
        plt.plot(time,x[:,2], 'blue', label='Recovered')
        plt.legend()
        plt.xlabel('Arbitrary time unit')
        plt.ylabel('Population numbers')
        string_for_title='SIROnGraph, k_SI ='+str(params['k_SI']) +', k_IR = '+ str(params['k_IR']) + ', S_ini = '+ str(params['S_ini']) + ', I_ini = '+str( params['I_ini']) + ', R_ini = '+ str(params['R_ini'])
        plt.title(string_for_title)
        string_for_savefig = 'SIROnGraph_graph_type_'+str(params['graph_type'])+'_k_SI='+str(params['k_SI']) +'_k_IR='+ str(params['k_IR']) + '_S_ini='+ str(params['S_ini']) + '_I_ini='+str( params['I_ini']) + '_R_ini='+ str(params['R_ini'])
        plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
        plt.close()
    
    if ('AllNodes' in params['plot_options']):
        print("Running: plot option: 'AllNodes'")
        ### --- loop over all saved data and plot SIRWellMixed on each region
        # use dic_x available in node attributes from graph to get 
        # number of nodes and number of reactions
        dic_x = nx.get_node_attributes(G,'x')        
        nbNodes = len(dic_x)
        nbReactions = len(dic_x[0])
        nodes = G.nodes()
        all_nodes_all_x_over_time = {}
        for n in nodes:
            x = np.zeros([nbTimeSteps, nbReactions])
            for i in range(0,nbTimeSteps):
                dic_x = outputs[i][1]            
                for j in range(0,nbReactions):
                    x[i,j] = dic_x[n][j]                    
            all_nodes_all_x_over_time[n] = x
            fig10 = plt.figure(figsize=(16, 10)) # in inches!
            plt.plot(time,x[:,0], 'red', label='Susceptible')
            plt.plot(time,x[:,1], 'green', label='Infected')
            plt.plot(time,x[:,2], 'blue', label='Recovered')
            plt.legend()
            plt.xlabel('Arbitrary time unit')
            plt.ylabel('Population numbers')
            string_for_title='SIROnGraph, region num = '+str(n)+', k_SI ='+str(params['k_SI']) +', k_IR = '+ str(params['k_IR']) + ', S_ini = '+ str(params['S_ini']) + ', I_ini = '+str( params['I_ini']) + ', R_ini = '+ str(params['R_ini'])
            plt.title(string_for_title)
            string_for_savefig = 'SIROnGraph_graph_type_'+str(params['graph_type'])+'_region_num='+str(n)+'_k_SI='+str(params['k_SI']) +'_k_IR='+ str(params['k_IR']) + '_S_ini='+ str(params['S_ini']) + '_I_ini='+str( params['I_ini']) + '_R_ini='+ str(params['R_ini'])
            plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
            plt.close()
                            
    if ('AllGraphs' in params['plot_options']):
        print("Running: plot option: 'AllGraphs'")            
        ### --- loop over all saved graph and plot at each saved time step    
        for j in range(0,nbTimeSteps):
            print("saving ",j," th time step = " ,time[j])        
            # dic_x = outputs[j][1]
            dic_x = copy.deepcopy(outputs[j][1])
            dic_positions = nx.get_node_attributes(G,'positions')        
            dic_S = {}
            dic_I = {}
            dic_R = {}
            for i in dic_x.keys():
                dic_S[i] = dic_x[i][0]
                dic_I[i] = dic_x[i][1]
                dic_R[i] = dic_x[i][2]
            
            list_S_values = [v for v in dic_S.values()]
            list_I_values = [v for v in dic_I.values()]
            list_R_values = [v for v in dic_R.values()]                        
            fig20 = plt.figure(figsize=(16, 10)) # in inches!  
            string_for_title='Susceptible population on Graph at time = '+str(round(time[j],6))
            plt.title(string_for_title)
            ### --- drawing nodes and edges separately so we can capture collection for colobar
            ec = nx.draw_networkx_edges(G, pos=dic_positions, alpha=0.2)
            nc = nx.draw_networkx_nodes(G, pos=dic_positions, nodelist=G.nodes(), node_color=list_S_values, with_labels=False, node_size=100 ,cmap=plt.cm.viridis)#, cmap=plt.cm.viridis)
            plt.colorbar(nc)
            plt.axis('off')    
            string_for_savefig = 'Susceptible_population_on_Graph_graph_type='+str(params['graph_type'])+'_'+str('%08d'%j)
            plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')        
            plt.close()
            
            fig30 = plt.figure(figsize=(16, 10)) # in inches!
            string_for_title='Infected population on Graph at time = '+str(round(time[j],6))
            plt.title(string_for_title)
            ### --- drawing nodes and edges separately so we can capture collection for colobar
            ec = nx.draw_networkx_edges(G, pos=dic_positions, alpha=0.2)
            nc = nx.draw_networkx_nodes(G, pos=dic_positions, nodelist=G.nodes(), node_color=list_I_values, with_labels=False, node_size=100 ,cmap=plt.cm.viridis)#, cmap=plt.cm.viridis)
            plt.colorbar(nc)
            plt.axis('off')    
            string_for_savefig = 'Infected_population_on_Graph_graph_type='+str(params['graph_type'])+'_'+str('%08d'%j)
            plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
            plt.close()
            
            fig40 = plt.figure(figsize=(16, 10)) # in inches!            
            string_for_title='Recovered population on Graph at time = '+str(round(time[j],6))
            plt.title(string_for_title)
            ### --- drawing nodes and edges separately so we can capture collection for colobar
            ec = nx.draw_networkx_edges(G, pos=dic_positions, alpha=0.2)
            nc = nx.draw_networkx_nodes(G, pos=dic_positions, nodelist=G.nodes(), node_color=list_R_values, with_labels=False, node_size=100 ,cmap=plt.cm.viridis)#, cmap=plt.cm.viridis)
            plt.colorbar(nc)
            plt.axis('off')
            string_for_savefig = 'Recovered_population_on_Graph_graph_type='+str(params['graph_type'])+'_'+str('%08d'%j)
            plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
            plt.close()
    if 'FinalGraph' in params['plot_options']:
        print("Running: plot option: 'FinalGraph'")
        # ### --- plot graph with final species update == state in Graph
        dic_positions = nx.get_node_attributes(G,'positions')        
        dic_x = nx.get_node_attributes(G,'x')
        dic_S = {}
        dic_I = {}
        dic_R = {}
        for i in dic_x.keys():
            dic_S[i] = dic_x[i][0]
            dic_I[i] = dic_x[i][1]
            dic_R[i] = dic_x[i][2]
        
        list_S_values = [v for v in dic_S.values()]
        list_I_values = [v for v in dic_I.values()]
        list_R_values = [v for v in dic_R.values()]
        fig20 = plt.figure(figsize=(16, 10)) # in inches!  
        string_for_title='Susceptible population on Graph at time = '+str(round(G.graph['time'],6))
        plt.title(string_for_title)
        # ### --- drawing nodes and edges separately so we can capture collection for colobar
        ec = nx.draw_networkx_edges(G, pos=dic_positions, alpha=0.2)
        nc = nx.draw_networkx_nodes(G, pos=dic_positions, nodelist=G.nodes(), node_color=list_S_values, with_labels=False, node_size=100 ,cmap=plt.cm.viridis)#, cmap=plt.cm.viridis)
        plt.colorbar(nc)
        plt.axis('off')    
        string_for_savefig = 'Susceptible_population_on_Final_Graph_at_time='+str(round(G.graph['time'],6))
        plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
        plt.close()
        
        fig30 = plt.figure(figsize=(16, 10)) # in inches!
        string_for_title='Infected population on Graph at time = '+str(round(G.graph['time'],6))
        plt.title(string_for_title)
        # ### --- drawing nodes and edges separately so we can capture collection for colobar
        ec = nx.draw_networkx_edges(G, pos=dic_positions, alpha=0.2)
        nc = nx.draw_networkx_nodes(G, pos=dic_positions, nodelist=G.nodes(), node_color=list_I_values, with_labels=False, node_size=100 ,cmap=plt.cm.viridis)#, cmap=plt.cm.viridis)
        plt.colorbar(nc)
        plt.axis('off')    
        string_for_savefig = 'Infected_population_on_Final_Graph_at_time='+str(round(G.graph['time'],6))
        plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
        plt.close()
        
        fig40 = plt.figure(figsize=(16, 10)) # in inches!            
        string_for_title='Recovered population on Graph at time = '+str(round(G.graph['time'],6))
        plt.title(string_for_title)
        # ### --- drawing nodes and edges separately so we can capture collection for colobar
        ec = nx.draw_networkx_edges(G, pos=dic_positions, alpha=0.2)
        nc = nx.draw_networkx_nodes(G, pos=dic_positions, nodelist=G.nodes(), node_color=list_R_values, with_labels=False, node_size=100 ,cmap=plt.cm.viridis)#, cmap=plt.cm.viridis)
        plt.colorbar(nc)
        plt.axis('off')
        string_for_savefig = 'Recovered_population_on_Final_Graph_at_time='+str(round(G.graph['time'],6))
        plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')        
        plt.close()

def plot_SIRDOnGraph(params,outputs, G):
    """
    plot plot_SIRDOnGraph model results
    S,I, R, D populations over time
    depending values contained in params['plot_options'] we can either:
    WellMixed: plot the well mixed limit by summing all compartments
    AllNodes: plot the dynamic on every node separately
    AllGraphs: plot all graphs over all time steps
    FinalGraph: plot only the final graph at the final time step
    """
    print('plot_SIRDOnGraph')
    # save outputs
    f = open(params['outputs_dir']+'/outputs.txt',"w")
    f.write(str(outputs))
    f.close()
    
    # save params inputs
    f = open(params['outputs_dir']+'/params.txt',"w")
    f.write(str(params))
    f.close()
    
    nbTimeSteps = len(outputs)
    time=np.zeros(nbTimeSteps)
    x = np.zeros([nbTimeSteps, params['nbStates']])
    dic_graphs = {}
    for i in range(0,nbTimeSteps):
        time[i] = outputs[i][0]
        dic_x = outputs[i][1]        
        for n in dic_x.keys():
            for j in range(0, params['nbStates']):
                x[i,j] += dic_x[n][j]
                                
    if ('WellMixed' in params['plot_options']):
        print("Running: plot option: 'WellMixed'")
        ### --- plot wellmixed model ( = sum over all compartments i.e. nodes of the graph)
        fig10 = plt.figure(figsize=(16, 10)) # in inches!
        plt.plot(time,x[:,0], 'red', label='Susceptible')
        plt.plot(time,x[:,1], 'green', label='Infected')
        plt.plot(time,x[:,2], 'blue', label='Recovered')
        plt.plot(time,x[:,3], 'black', label='Deaths')
        plt.legend()
        plt.xlabel('Arbitrary time unit')
        plt.ylabel('Population numbers')
        string_for_title='SIRDOnGraph, k_SI ='+str(params['k_SI']) +', k_IR = '+ str(params['k_IR'])+', k_ID = '+ str(params['k_ID']) + ', S_ini = '+ str(params['S_ini']) + ', I_ini = '+str( params['I_ini']) + ', R_ini = '+ str(params['R_ini'])+ ', D_ini = '+ str(params['D_ini'])
        plt.title(string_for_title)
        string_for_savefig = 'SIRDOnGraph_graph_type_'+str(params['graph_type'])+'_k_SI='+str(params['k_SI']) +'_k_IR='+ str(params['k_IR'])+'_k_ID='+ str(params['k_ID']) + '_S_ini='+ str(params['S_ini']) + '_I_ini='+str( params['I_ini']) + '_R_ini='+ str(params['R_ini'])+ '_D_ini='+ str(params['D_ini'])
        plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
        plt.close()
    
    if ('AllNodes' in params['plot_options']):
        print("Running: plot option: 'AllNodes'")
        ### --- loop over all saved data and plot SIRDWellMixed on each region
        # use dic_x available in node attributes from graph to get 
        # number of nodes and number of reactions
        dic_x = nx.get_node_attributes(G,'x')        
        nbNodes = len(dic_x)        
        nodes = G.nodes()
        all_nodes_all_x_over_time = {}
        for n in nodes:
            x = np.zeros([nbTimeSteps, params['nbStates']])
            for i in range(0,nbTimeSteps):
                dic_x = outputs[i][1]            
                for j in range(0,params['nbStates']):
                    x[i,j] = dic_x[n][j]  
                                      
            all_nodes_all_x_over_time[n] = x
            fig10 = plt.figure(figsize=(16, 10)) # in inches!
            plt.plot(time,x[:,0], 'red', label='Susceptible')
            plt.plot(time,x[:,1], 'green', label='Infected')
            plt.plot(time,x[:,2], 'blue', label='Recovered')
            plt.plot(time,x[:,3], 'black', label='Deaths')
            plt.legend()
            plt.xlabel('Arbitrary time unit')
            plt.ylabel('Population numbers')
            string_for_title='SIRDOnGraph, region num = '+str(n)+', k_SI ='+str(params['k_SI']) +', k_IR = '+ str(params['k_IR']) + ', S_ini = '+ str(params['S_ini']) + ', I_ini = '+str( params['I_ini']) + ', R_ini = '+ str(params['R_ini'])
            plt.title(string_for_title)
            string_for_savefig = 'SIRDOnGraph_graph_type_'+str(params['graph_type'])+'_region_num='+str(n)+'_k_SI='+str(params['k_SI']) +'_k_IR='+ str(params['k_IR'])+'_k_ID='+ str(params['k_ID']) + '_S_ini='+ str(params['S_ini']) + '_I_ini='+str( params['I_ini']) + '_R_ini='+ str(params['R_ini'])+ '_D_ini='+ str(params['D_ini'])
            plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
            plt.close()
                            
    if ('AllGraphs' in params['plot_options']):
        print("Running: plot option: 'AllGraphs'")            
        ### --- loop over all saved graph and plot at each saved time step    
        for j in range(0,nbTimeSteps):
            print("saving ",j," th time step = " ,time[j])        
            # dic_x = outputs[j][1]
            dic_x = copy.deepcopy(outputs[j][1])
            dic_positions = nx.get_node_attributes(G,'positions')        
            dic_S = {}
            dic_I = {}
            dic_R = {}
            dic_D = {}
            for i in dic_x.keys():
                dic_S[i] = dic_x[i][0]
                dic_I[i] = dic_x[i][1]
                dic_R[i] = dic_x[i][2]
                dic_D[i] = dic_x[i][3]
            
            list_S_values = [v for v in dic_S.values()]
            list_I_values = [v for v in dic_I.values()]
            list_R_values = [v for v in dic_R.values()]                        
            list_D_values = [v for v in dic_D.values()]                        
            fig20 = plt.figure(figsize=(16, 10)) # in inches!  
            string_for_title='Susceptible population on Graph at time = '+str(round(time[j],6))
            plt.title(string_for_title)
            ### --- drawing nodes and edges separately so we can capture collection for colobar
            ec = nx.draw_networkx_edges(G, pos=dic_positions, alpha=0.2)
            nc = nx.draw_networkx_nodes(G, pos=dic_positions, nodelist=G.nodes(), node_color=list_S_values, with_labels=False, node_size=100 ,cmap=plt.cm.viridis)#, cmap=plt.cm.viridis)
            plt.colorbar(nc)
            plt.axis('off')    
            string_for_savefig = 'Susceptible_population_on_Graph_graph_type='+str(params['graph_type'])+'_'+str('%08d'%j)
            plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')        
            plt.close()
            
            fig30 = plt.figure(figsize=(16, 10)) # in inches!
            string_for_title='Infected population on Graph at time = '+str(round(time[j],6))
            plt.title(string_for_title)
            ### --- drawing nodes and edges separately so we can capture collection for colobar
            ec = nx.draw_networkx_edges(G, pos=dic_positions, alpha=0.2)
            nc = nx.draw_networkx_nodes(G, pos=dic_positions, nodelist=G.nodes(), node_color=list_I_values, with_labels=False, node_size=100 ,cmap=plt.cm.viridis)#, cmap=plt.cm.viridis)
            plt.colorbar(nc)
            plt.axis('off')    
            string_for_savefig = 'Infected_population_on_Graph_graph_type='+str(params['graph_type'])+'_'+str('%08d'%j)
            plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
            plt.close()
            
            fig40 = plt.figure(figsize=(16, 10)) # in inches!            
            string_for_title='Recovered population on Graph at time = '+str(round(time[j],6))
            plt.title(string_for_title)
            ### --- drawing nodes and edges separately so we can capture collection for colobar
            ec = nx.draw_networkx_edges(G, pos=dic_positions, alpha=0.2)
            nc = nx.draw_networkx_nodes(G, pos=dic_positions, nodelist=G.nodes(), node_color=list_R_values, with_labels=False, node_size=100 ,cmap=plt.cm.viridis)#, cmap=plt.cm.viridis)
            plt.colorbar(nc)
            plt.axis('off')
            string_for_savefig = 'Recovered_population_on_Graph_graph_type='+str(params['graph_type'])+'_'+str('%08d'%j)
            plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
            plt.close()
            
            fig50 = plt.figure(figsize=(16, 10)) # in inches!            
            string_for_title='Deaths population on Graph at time = '+str(round(time[j],6))
            plt.title(string_for_title)
            ### --- drawing nodes and edges separately so we can capture collection for colobar
            ec = nx.draw_networkx_edges(G, pos=dic_positions, alpha=0.2)
            nc = nx.draw_networkx_nodes(G, pos=dic_positions, nodelist=G.nodes(), node_color=list_D_values, with_labels=False, node_size=100 ,cmap=plt.cm.viridis)#, cmap=plt.cm.viridis)
            plt.colorbar(nc)
            plt.axis('off')
            string_for_savefig = 'Deaths_population_on_Graph_graph_type='+str(params['graph_type'])+'_'+str('%08d'%j)
            plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
            plt.close()
            
    if 'FinalGraph' in params['plot_options']:
        print("Running: plot option: 'FinalGraph'")
        # ### --- plot graph with final species update == state in Graph
        dic_positions = nx.get_node_attributes(G,'positions')        
        dic_x = nx.get_node_attributes(G,'x')
        dic_S = {}
        dic_I = {}
        dic_R = {}
        dic_D = {}
        for i in dic_x.keys():
            dic_S[i] = dic_x[i][0]
            dic_I[i] = dic_x[i][1]
            dic_R[i] = dic_x[i][2]
            dic_D[i] = dic_x[i][3]
        
        list_S_values = [v for v in dic_S.values()]
        list_I_values = [v for v in dic_I.values()]
        list_R_values = [v for v in dic_R.values()]
        list_D_values = [v for v in dic_D.values()]
        fig20 = plt.figure(figsize=(16, 10)) # in inches!  
        string_for_title='Susceptible population on Graph at time = '+str(round(G.graph['time'],6))
        plt.title(string_for_title)
        # ### --- drawing nodes and edges separately so we can capture collection for colobar
        ec = nx.draw_networkx_edges(G, pos=dic_positions, alpha=0.2)
        nc = nx.draw_networkx_nodes(G, pos=dic_positions, nodelist=G.nodes(), node_color=list_S_values, with_labels=False, node_size=100 ,cmap=plt.cm.viridis)#, cmap=plt.cm.viridis)
        plt.colorbar(nc)
        plt.axis('off')    
        string_for_savefig = 'Susceptible_population_on_Final_Graph_at_time='+str(round(G.graph['time'],6))
        plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
        plt.close()
        
        fig30 = plt.figure(figsize=(16, 10)) # in inches!
        string_for_title='Infected population on Graph at time = '+str(round(G.graph['time'],6))
        plt.title(string_for_title)
        # ### --- drawing nodes and edges separately so we can capture collection for colobar
        ec = nx.draw_networkx_edges(G, pos=dic_positions, alpha=0.2)
        nc = nx.draw_networkx_nodes(G, pos=dic_positions, nodelist=G.nodes(), node_color=list_I_values, with_labels=False, node_size=100 ,cmap=plt.cm.viridis)#, cmap=plt.cm.viridis)
        plt.colorbar(nc)
        plt.axis('off')    
        string_for_savefig = 'Infected_population_on_Final_Graph_at_time='+str(round(G.graph['time'],6))
        plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
        plt.close()
        
        fig40 = plt.figure(figsize=(16, 10)) # in inches!            
        string_for_title='Recovered population on Graph at time = '+str(round(G.graph['time'],6))
        plt.title(string_for_title)
        # ### --- drawing nodes and edges separately so we can capture collection for colobar
        ec = nx.draw_networkx_edges(G, pos=dic_positions, alpha=0.2)
        nc = nx.draw_networkx_nodes(G, pos=dic_positions, nodelist=G.nodes(), node_color=list_R_values, with_labels=False, node_size=100 ,cmap=plt.cm.viridis)#, cmap=plt.cm.viridis)
        plt.colorbar(nc)
        plt.axis('off')
        string_for_savefig = 'Recovered_population_on_Final_Graph_at_time='+str(round(G.graph['time'],6))
        plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')        
        plt.close()
        
        fig50 = plt.figure(figsize=(16, 10)) # in inches!            
        string_for_title='Deaths population on Graph at time = '+str(round(G.graph['time'],6))
        plt.title(string_for_title)
        # ### --- drawing nodes and edges separately so we can capture collection for colobar
        ec = nx.draw_networkx_edges(G, pos=dic_positions, alpha=0.2)
        nc = nx.draw_networkx_nodes(G, pos=dic_positions, nodelist=G.nodes(), node_color=list_D_values, with_labels=False, node_size=100 ,cmap=plt.cm.viridis)#, cmap=plt.cm.viridis)
        plt.colorbar(nc)
        plt.axis('off')
        string_for_savefig = 'Deaths_population_on_Final_Graph_at_time='+str(round(G.graph['time'],6))
        plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')        
        plt.close()

def plot_SIRerlangOnGraph(params,outputs, G):
    """
    TODO
    """

def plot_TEST(params,outputs, G):
    """
    plot plot_SIRDOnGraph model results
    S,I, R, D populations over time
    depending values contained in params['plot_options'] we can either:
    WellMixed: plot the well mixed limit by summing all compartments
    AllNodes: plot the dynamic on every node separately
    AllGraphs: plot all graphs over all time steps
    FinalGraph: plot only the final graph at the final time step
    """
    print('plot_SIRDOnGraph')
    # save outputs
    f = open(params['outputs_dir']+'/outputs.txt',"w")
    f.write(str(outputs))
    f.close()
    
    # save params inputs
    f = open(params['outputs_dir']+'/params.txt',"w")
    f.write(str(params))
    f.close()
    
    nbTimeSteps = len(outputs)
    time=np.zeros(nbTimeSteps)
    x = np.zeros([nbTimeSteps, params['nbStates']])
    dic_graphs = {}
    for i in range(0,nbTimeSteps):
        time[i] = outputs[i][0]
        dic_x = outputs[i][1]        
        for n in dic_x.keys():
            for j in range(0, params['nbStates']):
                x[i,j] += dic_x[n][j]
                                
    if ('WellMixed' in params['plot_options']):
        print("Running: plot option: 'WellMixed'")
        ### --- plot wellmixed model ( = sum over all compartments i.e. nodes of the graph)
        fig10 = plt.figure(figsize=(16, 10)) # in inches!
        plt.plot(time,x[:,0], 'red', label='Susceptible')
        plt.plot(time,x[:,1], 'green', label='Infected')
        plt.plot(time,x[:,2], 'blue', label='Recovered')
        plt.plot(time,x[:,3], 'black', label='Deaths')
        plt.legend()
        plt.xlabel('Arbitrary time unit')
        plt.ylabel('Population numbers')
        string_for_title='TEST'
        plt.title(string_for_title)
        string_for_savefig = 'SIRDOnGraph_TEST_WellMixed'
        plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
        plt.close()
    
    if ('AllNodes' in params['plot_options']):
        print("Running: plot option: 'AllNodes'")
        ### --- loop over all saved data and plot SIRDWellMixed on each region
        # use dic_x available in node attributes from graph to get 
        # number of nodes and number of reactions
        dic_x = nx.get_node_attributes(G,'x')        
        nbNodes = len(dic_x)        
        nodes = G.nodes()
        all_nodes_all_x_over_time = {}
        for n in nodes:
            x = np.zeros([nbTimeSteps, params['nbStates']])
            for i in range(0,nbTimeSteps):
                dic_x = outputs[i][1]            
                for j in range(0,params['nbStates']):
                    x[i,j] = dic_x[n][j]  
                                      
            all_nodes_all_x_over_time[n] = x
            fig10 = plt.figure(figsize=(16, 10)) # in inches!
            plt.plot(time,x[:,0], 'red', label='Susceptible')
            plt.plot(time,x[:,1], 'green', label='Infected')
            plt.plot(time,x[:,2], 'blue', label='Recovered')
            plt.plot(time,x[:,3], 'black', label='Deaths')
            plt.legend()
            plt.xlabel('Arbitrary time unit')
            plt.ylabel('Population numbers')
            string_for_title='SIRDOnGraph_TEST'
            plt.title(string_for_title)
            string_for_savefig = 'SIRDOnGraph_TEST_node='+str(n)
            plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
            plt.close()
                            
    if ('AllGraphs' in params['plot_options']):
        print("Running: plot option: 'AllGraphs'")            
        ### --- loop over all saved graph and plot at each saved time step    
        for j in range(0,nbTimeSteps):
            print("saving ",j," th time step = " ,time[j])        
            # dic_x = outputs[j][1]
            dic_x = copy.deepcopy(outputs[j][1])
            dic_positions = nx.get_node_attributes(G,'positions')        
            dic_S = {}
            dic_I = {}
            dic_R = {}
            dic_D = {}
            for i in dic_x.keys():
                dic_S[i] = dic_x[i][0]
                dic_I[i] = dic_x[i][1]
                dic_R[i] = dic_x[i][2]
                dic_D[i] = dic_x[i][3]
            
            list_S_values = [v for v in dic_S.values()]
            list_I_values = [v for v in dic_I.values()]
            list_R_values = [v for v in dic_R.values()]                        
            list_D_values = [v for v in dic_D.values()]                        
            fig20 = plt.figure(figsize=(16, 10)) # in inches!  
            string_for_title='Susceptible population on Graph at time = '+str(round(time[j],6))
            plt.title(string_for_title)
            ### --- drawing nodes and edges separately so we can capture collection for colobar
            ec = nx.draw_networkx_edges(G, pos=dic_positions, alpha=0.2)
            nc = nx.draw_networkx_nodes(G, pos=dic_positions, nodelist=G.nodes(), node_color=list_S_values, with_labels=False, node_size=100 ,cmap=plt.cm.viridis)#, cmap=plt.cm.viridis)
            plt.colorbar(nc)
            plt.axis('off')    
            string_for_savefig = 'Susceptible_population_on_Graph_graph_type='+str(params['graph_type'])+'_'+str('%08d'%j)
            plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')        
            plt.close()
            
            fig30 = plt.figure(figsize=(16, 10)) # in inches!
            string_for_title='Infected population on Graph at time = '+str(round(time[j],6))
            plt.title(string_for_title)
            ### --- drawing nodes and edges separately so we can capture collection for colobar
            ec = nx.draw_networkx_edges(G, pos=dic_positions, alpha=0.2)
            nc = nx.draw_networkx_nodes(G, pos=dic_positions, nodelist=G.nodes(), node_color=list_I_values, with_labels=False, node_size=100 ,cmap=plt.cm.viridis)#, cmap=plt.cm.viridis)
            plt.colorbar(nc)
            plt.axis('off')    
            string_for_savefig = 'Infected_population_on_Graph_graph_type='+str(params['graph_type'])+'_'+str('%08d'%j)
            plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
            plt.close()
            
            fig40 = plt.figure(figsize=(16, 10)) # in inches!            
            string_for_title='Recovered population on Graph at time = '+str(round(time[j],6))
            plt.title(string_for_title)
            ### --- drawing nodes and edges separately so we can capture collection for colobar
            ec = nx.draw_networkx_edges(G, pos=dic_positions, alpha=0.2)
            nc = nx.draw_networkx_nodes(G, pos=dic_positions, nodelist=G.nodes(), node_color=list_R_values, with_labels=False, node_size=100 ,cmap=plt.cm.viridis)#, cmap=plt.cm.viridis)
            plt.colorbar(nc)
            plt.axis('off')
            string_for_savefig = 'Recovered_population_on_Graph_graph_type='+str(params['graph_type'])+'_'+str('%08d'%j)
            plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
            plt.close()
            
            fig50 = plt.figure(figsize=(16, 10)) # in inches!            
            string_for_title='Deaths population on Graph at time = '+str(round(time[j],6))
            plt.title(string_for_title)
            ### --- drawing nodes and edges separately so we can capture collection for colobar
            ec = nx.draw_networkx_edges(G, pos=dic_positions, alpha=0.2)
            nc = nx.draw_networkx_nodes(G, pos=dic_positions, nodelist=G.nodes(), node_color=list_D_values, with_labels=False, node_size=100 ,cmap=plt.cm.viridis)#, cmap=plt.cm.viridis)
            plt.colorbar(nc)
            plt.axis('off')
            string_for_savefig = 'Deaths_population_on_Graph_graph_type='+str(params['graph_type'])+'_'+str('%08d'%j)
            plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
            plt.close()
            
    if 'FinalGraph' in params['plot_options']:
        print("Running: plot option: 'FinalGraph'")
        # ### --- plot graph with final species update == state in Graph
        dic_positions = nx.get_node_attributes(G,'positions')        
        dic_x = nx.get_node_attributes(G,'x')
        dic_S = {}
        dic_I = {}
        dic_R = {}
        dic_D = {}
        for i in dic_x.keys():
            dic_S[i] = dic_x[i][0]
            dic_I[i] = dic_x[i][1]
            dic_R[i] = dic_x[i][2]
            dic_D[i] = dic_x[i][3]
        
        list_S_values = [v for v in dic_S.values()]
        list_I_values = [v for v in dic_I.values()]
        list_R_values = [v for v in dic_R.values()]
        list_D_values = [v for v in dic_D.values()]
        fig20 = plt.figure(figsize=(16, 10)) # in inches!  
        string_for_title='Susceptible population on Graph at time = '+str(round(G.graph['time'],6))
        plt.title(string_for_title)
        # ### --- drawing nodes and edges separately so we can capture collection for colobar
        ec = nx.draw_networkx_edges(G, pos=dic_positions, alpha=0.2)
        nc = nx.draw_networkx_nodes(G, pos=dic_positions, nodelist=G.nodes(), node_color=list_S_values, with_labels=False, node_size=100 ,cmap=plt.cm.viridis)#, cmap=plt.cm.viridis)
        plt.colorbar(nc)
        plt.axis('off')    
        string_for_savefig = 'Susceptible_population_on_Final_Graph_at_time='+str(round(G.graph['time'],6))
        plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
        plt.close()
        
        fig30 = plt.figure(figsize=(16, 10)) # in inches!
        string_for_title='Infected population on Graph at time = '+str(round(G.graph['time'],6))
        plt.title(string_for_title)
        # ### --- drawing nodes and edges separately so we can capture collection for colobar
        ec = nx.draw_networkx_edges(G, pos=dic_positions, alpha=0.2)
        nc = nx.draw_networkx_nodes(G, pos=dic_positions, nodelist=G.nodes(), node_color=list_I_values, with_labels=False, node_size=100 ,cmap=plt.cm.viridis)#, cmap=plt.cm.viridis)
        plt.colorbar(nc)
        plt.axis('off')    
        string_for_savefig = 'Infected_population_on_Final_Graph_at_time='+str(round(G.graph['time'],6))
        plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
        plt.close()
        
        fig40 = plt.figure(figsize=(16, 10)) # in inches!            
        string_for_title='Recovered population on Graph at time = '+str(round(G.graph['time'],6))
        plt.title(string_for_title)
        # ### --- drawing nodes and edges separately so we can capture collection for colobar
        ec = nx.draw_networkx_edges(G, pos=dic_positions, alpha=0.2)
        nc = nx.draw_networkx_nodes(G, pos=dic_positions, nodelist=G.nodes(), node_color=list_R_values, with_labels=False, node_size=100 ,cmap=plt.cm.viridis)#, cmap=plt.cm.viridis)
        plt.colorbar(nc)
        plt.axis('off')
        string_for_savefig = 'Recovered_population_on_Final_Graph_at_time='+str(round(G.graph['time'],6))
        plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')        
        plt.close()
        
        fig50 = plt.figure(figsize=(16, 10)) # in inches!            
        string_for_title='Deaths population on Graph at time = '+str(round(G.graph['time'],6))
        plt.title(string_for_title)
        # ### --- drawing nodes and edges separately so we can capture collection for colobar
        ec = nx.draw_networkx_edges(G, pos=dic_positions, alpha=0.2)
        nc = nx.draw_networkx_nodes(G, pos=dic_positions, nodelist=G.nodes(), node_color=list_D_values, with_labels=False, node_size=100 ,cmap=plt.cm.viridis)#, cmap=plt.cm.viridis)
        plt.colorbar(nc)
        plt.axis('off')
        string_for_savefig = 'Deaths_population_on_Final_Graph_at_time='+str(round(G.graph['time'],6))
        plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')        
        plt.close()

def plot_SIRWellMixedAgentBased(params,outputs):
    """
    plot SIRWellMixed model results with individual tracking of every agent
    S,I, R populations over time
    """
    # save outputs
    f = open(params['outputs_dir']+'/outputs.txt',"w")
    f.write(str(outputs))
    f.close()
    
    # save params inputs
    f = open(params['outputs_dir']+'/params.txt',"w")
    f.write(str(params))
    f.close()
    
    nbTimeSteps = len(outputs)
    time=np.zeros(nbTimeSteps)
    x = np.zeros([nbTimeSteps, 3])    
    for i in range(0,nbTimeSteps):
        time[i] = outputs[i][0]
        x[i,0] = len(outputs[i][1][0])
        x[i,1] = len(outputs[i][1][1])
        x[i,2] = len(outputs[i][1][2])
        
    fig2 = plt.figure(figsize=(16, 10)) # in inches!
    plt.plot(time,x[:,0], 'red', label='Susceptible')
    plt.plot(time,x[:,1], 'green', label='Infected')
    plt.plot(time,x[:,2], 'blue', label='Recovered')
    plt.legend()
    plt.xlabel('Arbitrary time unit')
    plt.ylabel('Population numbers')
    string_for_title='SIRWellMixedAgentBased, k_SI ='+str(params['k_SI']) +', k_IR = '+ str(params['k_IR']) + ', S_ini = '+ str(params['S_ini']) + ', I_ini = '+str( params['I_ini']) + ', R_ini = '+ str(params['R_ini'])
    plt.title(string_for_title)
    string_for_savefig = 'SIRWellMixedAgentBased_k_SI='+str(params['k_SI']) +'_k_IR='+ str(params['k_IR']) + '_S_ini='+ str(params['S_ini']) + '_I_ini='+str( params['I_ini']) + '_R_ini='+ str(params['R_ini'])
    plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
    # plt.close()
    plt.show()

def plot_SIerlangWellMixedAgentBased(params,outputs):
    """
    plot SIerlangWellMixed model results
    S,I
    respective indices for plot are 0, m, nbStates
    """
    # save outputs
    fname_outputs  = params['outputs_dir']+'/outputs_'+params['model_name']+'_m='+str(params['nb_intermediate_states'])+'.txt'
    f = open(fname_outputs,"w")
    f.write(str(outputs))
    f.close()
    
    # save params inputs
    fname_params  = params['outputs_dir']+'/params_'+params['model_name']+'_m='+str(params['nb_intermediate_states'])+'.txt'
    f = open(fname_params,"w")
    f.write(str(params))
    f.close()
    
    nbTimeSteps = len(outputs)
    time=np.zeros(nbTimeSteps)
    nbStates = params['nbStates']
    m = params["nb_intermediate_states"]
    x = np.zeros([nbTimeSteps, nbStates])
    fig2 = plt.figure(figsize=(16, 10)) # in inches!
    for i in range(0,nbTimeSteps):
        time[i] = outputs[i][0]
        for k in range(0,len(outputs[i][1])):
            x[i,k] = len(outputs[i][1][k])
    tmp_sum = np.zeros(nbTimeSteps)
    for ii in range(0,nbStates):
        if (ii ==0):
            _label='Susceptible'
            tmp_sum += x[:,ii]
        if (ii==m):
            _label='Infected'            
        if((ii>0) and (ii<m)):
            _label='I'+str(ii)     
            tmp_sum += x[:,ii]
        plt.plot(time,x[:,ii], label=_label)            
    
    plt.plot(time,tmp_sum,color='black',label='sum of S + I_i (i=1:m-1)')
    plt.legend()
    plt.xlabel('Arbitrary time unit')
    plt.ylabel('Population numbers')
    # plt.yscale('log')
    string_for_title='SIRerlangWellMixedAgentBased, k_SI ='+str(params['k_SI'])+' , nb_intermediate_states = ' +str(params['nb_intermediate_states']) + ', S_ini = '+ str(params['S_ini']) + ', I_ini = '+str( params['I_ini'])
    plt.title(string_for_title)
    string_for_savefig = 'SIRerlangWellMixedAgentBased_k_SI='+str(params['k_SI']) +'_nb_intermediate_states=' + str(params['nb_intermediate_states']) + '_S_ini='+ str(params['S_ini']) + '_I_ini='+str( params['I_ini'])
    plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
    # plt.close() 
    plt.show()  
    # plot waiting times    
    # plot waiting times    
    if ('WaitingTimes' in params['plot_options']):
        I_m = x[:,nbStates-1]        
        x_ini = params['x_ini']
        S_ini = x_ini[0]
        I_m_ini   = x_ini[nbStates-1]    
        
        # recover arrival times from L to I
        arrival_times_from_S_to_I_1 = []
        arrival_times_from_S_to_I_m = []
        for k in range(0,len(S_ini)):
            s = S_ini[k] # get id
            found_arrival_I_1 = False
            found_arrival_I_m = False
            for i in range(0,nbTimeSteps):                
                I_1_list = outputs[i][1][1]
                I_m_list = outputs[i][1][nbStates-1]
                _time = outputs[i][0]
                if (s in I_1_list) and not(found_arrival_I_1):
                    arrival_times_from_S_to_I_1.append(_time)
                    found_arrival_I_1 = True
                if (s in I_m_list) and not(found_arrival_I_m):
                    arrival_times_from_S_to_I_m.append(_time)
                    found_arrival_I_m = True
        
        arrival_times_from_I_1_to_I_m_TEST = []
        for kk in range(0,len(arrival_times_from_S_to_I_m)):
            arrival_times_from_I_1_to_I_m_TEST.append(arrival_times_from_S_to_I_m[kk] - arrival_times_from_S_to_I_1[kk])
        
        # recover arrival times from L to I
        arrival_times_from_I_1_to_I_m = []
        found_ids_in_I_1 = []
        for i in range(0,nbTimeSteps):
            I_1_list = outputs[i][1][1]            
            for I_1_id in I_1_list:
                if not(I_1_id in found_ids_in_I_1):
                    found_arrival = False
                    time_start = outputs[i][0]
                    for j in range(i,nbTimeSteps):
                        I_m_list = outputs[j][1][nbStates-1]
                        if (I_1_id in I_m_list) and not(found_arrival):
                            time_end = outputs[j][0]
                            arrival_times_from_I_1_to_I_m.append(time_end-time_start)
                            found_arrival = True
                            found_ids_in_I_1.append(I_1_id)
        # convert to numpy        
        arrival_times_from_S_to_I_1   =  np.array(arrival_times_from_S_to_I_1)
        arrival_times_from_I_1_to_I_m =  np.array(arrival_times_from_I_1_to_I_m)
        # arrival_times_from_I_1_to_I_m_TEST = np.array(arrival_times_from_I_1_to_I_m_TEST)
        # print('arrival_times_from_I_1_to_I_m_TEST = ' , arrival_times_from_I_1_to_I_m_TEST)
        print('arrival_times_from_I_1_to_I_m = ' , arrival_times_from_I_1_to_I_m)
        
        
        # m * rate gives rate of erlang distrib
        k_SI = params['k_SI']
        # k_LI = params['k_LI']
        m = params['nb_intermediate_states']
        m_tmp = m - 1
        Lambda_SI_tmp = m_tmp * k_SI * len(S_ini)
        time_real_I_1_to_I_m = np.linspace(0, np.max(arrival_times_from_I_1_to_I_m), 1000)
        time_real_S_to_I_m = np.linspace(0, np.max(arrival_times_from_S_to_I_m), 1000)
        # compute (m-1)!
        
        facto_m_tmp_minus_1 = 1
        for i in range(1,m_tmp):
            facto_m_tmp_minus_1 *= i
        erlang_distrib_I_1_to_I_m = Lambda_SI_tmp * (Lambda_SI_tmp * time_real_I_1_to_I_m) **(m_tmp-1) * np.exp(-Lambda_SI_tmp * time_real_I_1_to_I_m) / facto_m_tmp_minus_1        
        
        facto_m_minus_1 = 1
        for i in range(1,m):
            facto_m_minus_1 *= i
        Lambda_SI = m * k_SI * len(S_ini)
        erlang_distrib_S_to_I_m = Lambda_SI * (Lambda_SI * time_real_S_to_I_m) **(m-1) * np.exp(-Lambda_SI * time_real_S_to_I_m) / facto_m_minus_1        
        
        # plot histogram for waiting time distribution from L to I_m
        fig100 = plt.figure(figsize=(16, 10)) # in inches!
        plt.hist(arrival_times_from_I_1_to_I_m, bins=int(10), density=True, label='Waiting times distribution (from I_1 to I_m)')
        plt.plot(time_real_I_1_to_I_m, erlang_distrib_I_1_to_I_m, label='Erlang distribution (from I_1 to I_m)')
        plt.xlabel('Waiting times')
        plt.ylabel('PDF')                
        plt.legend()
        string_for_title='SIAndIntermediateErlangWellMixedAgentBased Waiting times from I_1 to I_m, nb_intermediate_states = ' +str(params['nb_intermediate_states']) + ', S_ini = '+ str(params['S_ini']) +  ', I_ini = '+str( params['I_ini'])
        plt.title(string_for_title)
        string_for_savefig = 'SIAndIntermediateErlangWellMixedAgentBased_waiting_times_pdf_from_I_1_to_I_m_k_SI='+str(params['k_SI'])+'_nb_intermediate_states=' + str(params['nb_intermediate_states']) + '_S_ini='+ str(params['S_ini'])+ '_I_ini='+str( params['I_ini'])
        plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
        # plt.close()
        plt.show()
        
        # plot histogram for waiting time distribution from S to L
        fig100 = plt.figure(figsize=(16, 10)) # in inches!
        plt.hist(arrival_times_from_S_to_I_1, bins=int(10), density=True, label='Waiting times distribution (from S to I_1)')
        plt.xlabel('Waiting times')
        plt.ylabel('PDF')                
        plt.legend()
        string_for_title='SIAndIntermediateErlangWellMixedAgentBased Waiting times from S to I_1, nb_intermediate_states = ' +str(params['nb_intermediate_states']) + ', S_ini = '+ str(params['S_ini']) + ', I_ini = '+str( params['I_ini'])
        plt.title(string_for_title)
        string_for_savefig = 'SIAndIntermediateErlangWellMixedAgentBased_waiting_times_pdf_from_S_to_I_1_k_SI='+str(params['k_SI']) +'_nb_intermediate_states=' + str(params['nb_intermediate_states']) + '_S_ini='+ str(params['S_ini']) + '_I_ini='+str( params['I_ini'])
        plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
        # plt.close()
        plt.show()
        
        # plot histogram for waiting time distribution from S to L
        fig100 = plt.figure(figsize=(16, 10)) # in inches!
        plt.hist(arrival_times_from_S_to_I_m, bins=int(10), density=True, label='Waiting times distribution (from S to I_m)')
        plt.plot(time_real_S_to_I_m, erlang_distrib_S_to_I_m, label='Erlang distribution (from S to I_m)')
        plt.xlabel('Waiting times')
        plt.ylabel('PDF')                
        plt.legend()
        string_for_title='SIAndIntermediateErlangWellMixedAgentBased Waiting times from S to I_m, nb_intermediate_states = ' +str(params['nb_intermediate_states']) + ', S_ini = '+ str(params['S_ini']) + ', I_ini = '+str( params['I_ini'])
        plt.title(string_for_title)
        string_for_savefig = 'SIAndIntermediateErlangWellMixedAgentBased_waiting_times_pdf_from_S_to_I_m_k_SI='+str(params['k_SI']) +'_nb_intermediate_states=' + str(params['nb_intermediate_states']) + '_S_ini='+ str(params['S_ini'])+ '_I_ini='+str( params['I_ini'])
        plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
        # plt.close()
        plt.show()


def plot_SIerlangModifiedWellMixedAgentBased(params,outputs, waiting_times):
    """
    plot SIerlangModifiedWellMixedAgentBased model results
    S,I
    respective indices for plot are 0, m, nbStates
    """
    # save outputs
    fname_outputs  = params['outputs_dir']+'/outputs_'+params['model_name']+'_m='+str(params['nb_intermediate_states'])+'.txt'
    f = open(fname_outputs,"w")
    f.write(str(outputs))
    f.close()
    
    # save params inputs
    fname_params  = params['outputs_dir']+'/params_'+params['model_name']+'_m='+str(params['nb_intermediate_states'])+'.txt'
    f = open(fname_params,"w")
    f.write(str(params))
    f.close()
    
    nbTimeSteps = len(outputs)
    time=np.zeros(nbTimeSteps)
    nbStates = params['nbStates']
    m = params["nb_intermediate_states"]
    x = np.zeros([nbTimeSteps, nbStates])
    fig2 = plt.figure(figsize=(16, 10)) # in inches!
    for i in range(0,nbTimeSteps):
        time[i] = outputs[i][0]
        for k in range(0,len(outputs[i][1])):
            x[i,k] = len(outputs[i][1][k])
    tmp_sum = np.zeros(nbTimeSteps)
    for ii in range(0,nbStates):
        if (ii ==0):
            _label='Susceptible'
            tmp_sum += x[:,ii]
        if (ii==m):
            _label='Infected'            
        if((ii>0) and (ii<m)):
            _label='I'+str(ii)     
            tmp_sum += x[:,ii]
        plt.plot(time,x[:,ii], label=_label)            
    
    plt.plot(time,tmp_sum,color='black',label='sum of S + I_i (i=1:m-1)')
    plt.legend()
    plt.xlabel('Arbitrary time unit')
    plt.ylabel('Population numbers')
    # plt.yscale('log')
    string_for_title='SIerlangModifiedWellMixedAgentBased, k_SI ='+str(params['k_SI'])+' , nb_intermediate_states = ' +str(params['nb_intermediate_states']) + ', S_ini = '+ str(params['S_ini']) + ', I_ini = '+str( params['I_ini'])
    plt.title(string_for_title)
    string_for_savefig = 'SIerlangModifiedWellMixedAgentBased_k_SI='+str(params['k_SI']) +'_nb_intermediate_states=' + str(params['nb_intermediate_states']) + '_S_ini='+ str(params['S_ini']) + '_I_ini='+str( params['I_ini'])
    plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
    # plt.close() 
    plt.show()  
    # plot waiting times    
    
    if ('WaitingTimes' in params['plot_options']):
        I_m = x[:,nbStates-1]        
        x_ini = params['x_ini']
        S_ini = x_ini[0]
        I_m_ini   = x_ini[nbStates-1]    
        
        arrival_times_from_S_to_I_m = []        
        arrival_times_from_S_to_I_m = waiting_times
        arrival_times_from_S_to_I_m   =  np.array(arrival_times_from_S_to_I_m)
        
        # m * rate gives rate of erlang distrib
        k_SI = params['k_SI']        
        m = params['nb_intermediate_states']        
        time_real_S_to_I_m = np.linspace(0, np.max(arrival_times_from_S_to_I_m), 1000)
        
        # compute (m-1)!
        facto_m_minus_1 = 1
        for i in range(1,m):
            facto_m_minus_1 *= i
            
        Lambda_SI = m * k_SI
        erlang_distrib_S_to_I_m = Lambda_SI * (Lambda_SI * time_real_S_to_I_m) **(m-1) * np.exp(-Lambda_SI * time_real_S_to_I_m) / facto_m_minus_1        
        
        
        # plot histogram for waiting time distribution from S to I_m
        fig100 = plt.figure(figsize=(16, 10)) # in inches!
        plt.hist(arrival_times_from_S_to_I_m, bins=int(30), density=True, label='Waiting times distribution (from S to I_m)')
        plt.plot(time_real_S_to_I_m, erlang_distrib_S_to_I_m, label='Erlang distribution (from S to I_m)')
        plt.xlabel('Waiting times')
        plt.ylabel('PDF')                
        plt.legend()
        string_for_title='SIAndIntermediateErlangWellMixedAgentBased Waiting times from S to I_m, nb_intermediate_states = ' +str(params['nb_intermediate_states']) + ', S_ini = '+ str(params['S_ini']) + ', I_ini = '+str( params['I_ini'])
        plt.title(string_for_title)
        string_for_savefig = 'SIAndIntermediateErlangWellMixedAgentBased_waiting_times_pdf_from_S_to_I_m_k_SI='+str(params['k_SI']) +'_nb_intermediate_states=' + str(params['nb_intermediate_states']) + '_S_ini='+ str(params['S_ini'])+ '_I_ini='+str( params['I_ini'])
        plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
        # plt.close()
        plt.show()


def plot_SIAndIntermediateErlangWellMixedAgentBased(params,outputs):
    """
    plot SIerlangWellMixed model results
    S,I
    respective indices for plot are 0, m, nbStates
    """
        # save outputs
    fname_outputs  = params['outputs_dir']+'/outputs_'+params['model_name']+'_m='+str(params['nb_intermediate_states'])+'.txt'
    f = open(fname_outputs,"w")
    f.write(str(outputs))
    f.close()
    
    # save params inputs
    fname_params  = params['outputs_dir']+'/params_'+params['model_name']+'_m='+str(params['nb_intermediate_states'])+'.txt'
    f = open(fname_params,"w")
    f.write(str(params))
    f.close()
    
    nbTimeSteps = len(outputs)
    time=np.zeros(nbTimeSteps)
    nbStates = params['nbStates']
    m = params["nb_intermediate_states"]
    x = np.zeros([nbTimeSteps, nbStates])
    fig2 = plt.figure(figsize=(16, 10)) # in inches!
    for i in range(0,nbTimeSteps):
        time[i] = outputs[i][0]        
        for k in range(0,len(outputs[i][1])):
            x[i,k] = len(outputs[i][1][k])
    tmp_sum = np.zeros(nbTimeSteps)
    for ii in range(0,nbStates):
        if (ii ==0):
            _label='Susceptible'
            tmp_sum += x[:,ii]
        if (ii ==1):
            _label='Intermediate L'
            tmp_sum += x[:,ii]
        if (ii==params['nbStates']-1):
            _label='Infected'     
        if((ii>1) and (ii<params['nbStates']-1)):
            _label='I'+str(ii-1)     
            tmp_sum += x[:,ii]
        plt.plot(time,x[:,ii], label=_label)            
    
    plt.plot(time,tmp_sum,color='black',label='sum of S + L + I_i (i=1:m-1)')
    plt.legend()
    plt.xlabel('Arbitrary time unit')
    plt.ylabel('Population numbers')
    string_for_title='SIAndIntermediateErlangWellMixedAgentBased, k_SL ='+str(params['k_SL'])+' k_LI ='+str(params['k_LI'])+' , nb_intermediate_states = ' +str(params['nb_intermediate_states']) + ', S_ini = '+ str(params['S_ini']) + ', I_ini = '+str( params['I_ini'])
    plt.title(string_for_title)
    string_for_savefig = 'SIAndIntermediateErlangWellMixedAgentBased_k_SL='+str(params['k_SL'])+str(params['k_LI'])+'_k_LI=' +'_nb_intermediate_states=' + str(params['nb_intermediate_states']) + '_S_ini='+ str(params['S_ini']) + '_I_ini='+str( params['I_ini'])
    plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
    # plt.close()
    plt.show()
    
    # plot waiting times    
    if ('WaitingTimes' in params['plot_options']):
        I_m = x[:,nbStates-1]        
        x_ini = params['x_ini']
        S_ini = x_ini[0]
        I_m_ini   = x_ini[nbStates-1]    
        
        # recover arrival times from L to I
        arrival_times_from_S_to_L = []
        for k in range(0,len(S_ini)):
            s = S_ini[k] # get id
            found_arrival = False
            for i in range(0,nbTimeSteps):                
                L_list = outputs[i][1][1]
                _time = outputs[i][0]
                if (s in L_list) and not(found_arrival):
                    arrival_times_from_S_to_L.append(_time)
                    found_arrival = True
        
        # recover arrival times from L to I
        arrival_times_from_L_to_I_m = []
        found_ids_in_L = []
        for i in range(0,nbTimeSteps):
            L_list = outputs[i][1][1]            
            for L_id in L_list:
                if not(L_id in found_ids_in_L):
                    found_arrival = False
                    time_start = outputs[i][0]
                    for j in range(i,nbTimeSteps):
                        I_m_list = outputs[j][1][nbStates-1]
                        if (L_id in I_m_list) and not(found_arrival):
                            time_end = outputs[j][0]
                            arrival_times_from_L_to_I_m.append(time_end-time_start)
                            found_arrival = True
                            found_ids_in_L.append(L_id)
        # convert to numpy        
        arrival_times_from_S_to_L   =  np.array(arrival_times_from_S_to_L)
        arrival_times_from_L_to_I_m =  np.array(arrival_times_from_L_to_I_m)
            
        # m * rate gives rate of erlang distrib
        k_SL = params['k_SL']
        k_LI = params['k_LI']
        m = params['nb_intermediate_states']
        Lambda_LI = m * k_LI
        time_real_L_to_I = np.linspace(0, np.max(arrival_times_from_L_to_I_m), 1000)
        time_real_S_to_L = np.linspace(0, np.max(arrival_times_from_S_to_L), 1000)
        # compute (m-1)!
        facto_m_minus_1 = 1
        for i in range(1,m):
            facto_m_minus_1 *= i
        erlang_distrib_L_to_I_m = Lambda_LI * (Lambda_LI * time_real_L_to_I) **(m-1) * np.exp(-Lambda_LI * time_real_L_to_I) / facto_m_minus_1        
        
        # expo distribution from S to L
        m_tmp=1
        Lambda_SL = k_SL
        erlang_distrib_S_to_L = Lambda_SL * (Lambda_SL * time_real_S_to_L) **(m_tmp-1) * np.exp(-Lambda_SL * time_real_S_to_L) / 1.   
                
        # plot histogram for waiting time distribution from L to I_m
        fig100 = plt.figure(figsize=(16, 10)) # in inches!
        plt.hist(arrival_times_from_L_to_I_m, bins=int(10), density=True, label='Waiting times distribution (from L to I_m)')
        plt.plot(time_real_L_to_I, erlang_distrib_L_to_I_m, label='Erlang distribution (from L to I_m)')
        plt.xlabel('Waiting times')
        plt.ylabel('PDF')                
        plt.legend()
        string_for_title='SIAndIntermediateErlangWellMixedAgentBased Waiting times from L to I, nb_intermediate_states = ' +str(params['nb_intermediate_states']) + ', S_ini = '+ str(params['S_ini']) + ', L_ini = '+str( params['L_ini'])+ ', I_ini = '+str( params['I_ini'])
        plt.title(string_for_title)
        string_for_savefig = 'SIAndIntermediateErlangWellMixedAgentBased_waiting_times_pdf_from_L_to_I_k_SL='+str(params['k_SL']) +'_k_LI='+str(params['k_LI']) +'_nb_intermediate_states=' + str(params['nb_intermediate_states']) + '_S_ini='+ str(params['S_ini'])+ '_L_ini='+str( params['L_ini']) + '_I_ini='+str( params['I_ini'])
        plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
        # plt.close()
        plt.show()
        
        # plot histogram for waiting time distribution from S to L
        fig100 = plt.figure(figsize=(16, 10)) # in inches!
        plt.hist(arrival_times_from_S_to_L, bins=int(10), density=True, label='Waiting times distribution (from S to L)')
        plt.plot(time_real_S_to_L, erlang_distrib_S_to_L, label='Erlang distribution (from S to L)')
        plt.xlabel('Waiting times')
        plt.ylabel('PDF')                
        plt.legend()
        string_for_title='SIAndIntermediateErlangWellMixedAgentBased Waiting times from S to L, nb_intermediate_states = ' +str(params['nb_intermediate_states']) + ', S_ini = '+ str(params['S_ini']) + ', L_ini = '+str( params['L_ini'])+ ', I_ini = '+str( params['I_ini'])
        plt.title(string_for_title)
        string_for_savefig = 'SIAndIntermediateErlangWellMixedAgentBased_waiting_times_pdf_from_S_to_L_k_SL='+str(params['k_SL']) +'_k_LI='+str(params['k_LI']) +'_nb_intermediate_states=' + str(params['nb_intermediate_states']) + '_S_ini='+ str(params['S_ini'])+ '_L_ini='+str( params['L_ini']) + '_I_ini='+str( params['I_ini'])
        plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
        # plt.close()
        plt.show()


def plot_SIRerlangKrylovaWellMixedAgentBased(params,outputs):
    """
    plot SIRerlangKrylovaWellMixedAgentBased model results
    S,I_1,...I_m, R
    respective indices for plot are 0, m, nbStates
    """
    # save outputs
    fname_outputs  = params['outputs_dir']+'/outputs_'+params['model_name']+'_m='+str(params['nb_intermediate_states'])+'.txt'
    f = open(fname_outputs,"w")
    f.write(str(outputs))
    f.close()
    
    # save params inputs
    fname_params  = params['outputs_dir']+'/params_'+params['model_name']+'_m='+str(params['nb_intermediate_states'])+'.txt'
    f = open(fname_params,"w")
    f.write(str(params))
    f.close()
    
    nbTimeSteps = len(outputs)
    time=np.zeros(nbTimeSteps)
    nbStates = params['nbStates']
    nbReactions = params['nbReactions']
    m = params["nb_intermediate_states"]
    x = np.zeros([nbTimeSteps, nbStates])
    fig2 = plt.figure(figsize=(16, 10)) # in inches!
    for i in range(0,nbTimeSteps):
        time[i] = outputs[i][0]
        for k in range(0,len(outputs[i][1])):
            x[i,k] = len(outputs[i][1][k])
    tmp_sum = np.zeros(nbTimeSteps)
    for ii in range(0,nbStates):
        if (ii ==0):
            _label='Susceptible'
            tmp_sum += x[:,ii]
        if((ii>0) and (ii<=nbStates-2)):
            _label='I'+str(ii)     
            tmp_sum += x[:,ii]
        if (ii == nbStates -1):
            _label ='R'
        plt.plot(time,x[:,ii], label=_label)            
    
    plt.plot(time,tmp_sum,color='black',label='sum of S + I_i (i=1:m-1)')
    plt.legend()
    plt.xlabel('Arbitrary time unit')
    plt.ylabel('Population numbers')
    # plt.yscale('log')
    string_for_title='SIRerlangKrylovaWellMixedAgentBased, k_SI ='+str(params['k_SI'])+', k_IR ='+str(params['k_IR'])+' , nb_intermediate_states = ' +str(params['nb_intermediate_states']) + ', S_ini = '+ str(params['S_ini']) + ', I_ini = '+str( params['I_ini'])
    plt.title(string_for_title)
    string_for_savefig = 'SIRerlangKrylovaWellMixedAgentBased_k_SI='+str(params['k_SI']) +'_k_IR='+str(params['k_IR'])+'_nb_intermediate_states=' + str(params['nb_intermediate_states']) + '_S_ini='+ str(params['S_ini']) + '_I_ini='+str( params['I_ini'])
    plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
    # plt.close() 
    plt.show()  
    # plot waiting times    
    # plot waiting times    
    if ('WaitingTimes' in params['plot_options']):
        I_m = x[:,nbStates-1]        
        x_ini = params['x_ini']
        S_ini = x_ini[0]
        I_m_ini   = x_ini[nbStates-1]    
        
        # # recover arrival times from L to I
        # arrival_times_from_S_to_I_1 = []
        # arrival_times_from_S_to_I_m = []
        # for k in range(0,len(S_ini)):
            # s = S_ini[k] # get id
            # found_arrival_I_1 = False
            # found_arrival_I_m = False
            # for i in range(0,nbTimeSteps):                
                # I_1_list = outputs[i][1][1]
                # I_m_list = outputs[i][1][nbStates-1]
                # _time = outputs[i][0]
                # if (s in I_1_list) and not(found_arrival_I_1):
                    # arrival_times_from_S_to_I_1.append(_time)
                    # found_arrival_I_1 = True
                # if (s in I_m_list) and not(found_arrival_I_m):
                    # arrival_times_from_S_to_I_m.append(_time)
                    # found_arrival_I_m = True
        
        # arrival_times_from_I_1_to_I_m_TEST = []
        # for kk in range(0,len(arrival_times_from_S_to_I_m)):
            # arrival_times_from_I_1_to_I_m_TEST.append(arrival_times_from_S_to_I_m[kk] - arrival_times_from_S_to_I_1[kk])
        
        # recover arrival times from L to I
        arrival_times_from_I_1_to_R = []
        found_ids_in_I_1 = []
        for i in range(0,nbTimeSteps):
            I_1_list = outputs[i][1][1]            
            for I_1_id in I_1_list:
                if not(I_1_id in found_ids_in_I_1):
                    found_arrival = False
                    time_start = outputs[i][0]
                    for j in range(i,nbTimeSteps):
                        R_list = outputs[j][1][nbStates-1]
                        if (I_1_id in R_list) and not(found_arrival):
                            time_end = outputs[j][0]
                            arrival_times_from_I_1_to_R.append(time_end-time_start)
                            found_arrival = True
                            found_ids_in_I_1.append(I_1_id)
        # convert to numpy        
        # arrival_times_from_S_to_I_1   =  np.array(arrival_times_from_S_to_I_1)
        arrival_times_from_I_1_to_R =  np.array(arrival_times_from_I_1_to_R)
        
        # m * rate gives rate of erlang distrib
        k_SI = params['k_SI']
        k_IR = params['k_IR']
        m = params['nb_intermediate_states']        
        time_real_I_1_to_R = np.linspace(0, np.max(arrival_times_from_I_1_to_R), 1000)        
        # compute (m-1)!
        
        facto_m_minus_1 = 1
        for i in range(1,m):
            facto_m_minus_1 *= i
        Lamba_IR = m * k_IR# * len(S_ini)
        erlang_distrib_I_1_to_R = Lamba_IR * (Lamba_IR * time_real_I_1_to_R) **(m-1) * np.exp(-Lamba_IR * time_real_I_1_to_R) / facto_m_minus_1        
        
        # plot histogram for waiting time distribution from I_1 to R
        fig100 = plt.figure(figsize=(16, 10)) # in inches!
        plt.hist(arrival_times_from_I_1_to_R, bins=int(50), density=True, label='Waiting times distribution (from I_1 to R)')
        plt.plot(time_real_I_1_to_R, erlang_distrib_I_1_to_R, label='Erlang distribution (from I_1 to R)')
        plt.xlabel('Waiting times')
        plt.ylabel('PDF')                
        plt.legend()
        string_for_title='SIRerlangKrylovaWellMixedAgentBased Waiting times from I_1 to R, nb_intermediate_states = ' +str(params['nb_intermediate_states']) + ', S_ini = '+ str(params['S_ini']) +  ', I_ini = '+str( params['I_ini'])
        plt.title(string_for_title)
        string_for_savefig = 'SIRerlangKrylovaWellMixedAgentBased_times_pdf_from_I_1_to_R_k_SI='+str(params['k_SI'])+'_k_IR='+str(params['k_IR'])+'_nb_intermediate_states=' + str(params['nb_intermediate_states']) + '_S_ini='+ str(params['S_ini'])+ '_I_ini='+str( params['I_ini'])
        plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
        # plt.close()
        plt.show()
        
def plot_SIROnGraphAgentBased(params,outputs, G):
    """
    plot plot_SIROnGraph model results
    S,I, R populations over time
    depending values contained in params['plot_options'] we can either:
    WellMixed: plot the well mixed limit by summing all compartments
    AllNodes: plot the dynamic on every node separately
    AllGraphs: plot all graphs over all time steps
    FinalGraph: plot only the final graph at the final time step
    """
    print('plot_SIROnGraphAgentBased')
    # save outputs
    f = open(params['outputs_dir']+'/outputs.txt',"w")
    f.write(str(outputs))
    f.close()
    
    # save params inputs
    f = open(params['outputs_dir']+'/params.txt',"w")
    f.write(str(params))
    f.close()
    
    nbTimeSteps = len(outputs)
    mytime=np.zeros(nbTimeSteps)
    x = np.zeros([nbTimeSteps, 3])
    dic_graphs = {}
    for i in range(0,nbTimeSteps):
        mytime[i] = outputs[i][0]
        dic_x = outputs[i][1]        
        for n in dic_x.keys():
            x[i,0] += len(dic_x[n][0])
            x[i,1] += len(dic_x[n][1])
            x[i,2] += len(dic_x[n][2])      
    
    if ('WellMixed' in params['plot_options']):
        print("Running: plot option: 'WellMixed'")
        ### --- plot wellmixed model ( = sum over all compartments i.e. nodes of the graph)
        fig10 = plt.figure(figsize=(16, 10)) # in inches!
        plt.plot(mytime,x[:,0], 'red', label='Susceptible')
        plt.plot(mytime,x[:,1], 'green', label='Infected')
        plt.plot(mytime,x[:,2], 'blue', label='Recovered')
        plt.legend()
        plt.xlabel('Arbitrary time unit')
        plt.ylabel('Population numbers')
        string_for_title='SIROnGraphAgentBased, k_SI ='+str(params['k_SI']) +', k_IR = '+ str(params['k_IR']) + ', S_ini = '+ str(params['S_ini']) + ', I_ini = '+str( params['I_ini']) + ', R_ini = '+ str(params['R_ini'])
        plt.title(string_for_title)
        string_for_savefig = 'SIROnGraphAgentBased_graph_type_'+str(params['graph_type'])+'_k_SI='+str(params['k_SI']) +'_k_IR='+ str(params['k_IR']) + '_S_ini='+ str(params['S_ini']) + '_I_ini='+str( params['I_ini']) + '_R_ini='+ str(params['R_ini'])
        plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
        plt.close()
    
    if ('AllNodes' in params['plot_options']):
        print("Running: plot option: 'AllNodes'")
        ### --- loop over all saved data and plot SIRWellMixed on each region
        # use dic_x available in node attributes from graph to get 
        # number of nodes and number of reactions
        dic_x = nx.get_node_attributes(G,'x')        
        nbNodes = len(dic_x)
        nbStates = len(dic_x[0])
        nodes = G.nodes()
        all_nodes_all_x_over_time = {}
        for n in nodes:
            x = np.zeros([nbTimeSteps, nbStates])
            for i in range(0,nbTimeSteps):
                dic_x = outputs[i][1]            
                for j in range(0,nbStates):
                    x[i,j] = len(dic_x[n][j])
            all_nodes_all_x_over_time[n] = x
            fig10 = plt.figure(figsize=(16, 10)) # in inches!
            plt.plot(mytime,x[:,0], 'red', label='Susceptible')
            plt.plot(mytime,x[:,1], 'green', label='Infected')
            plt.plot(mytime,x[:,2], 'blue', label='Recovered')
            plt.legend()
            plt.xlabel('Arbitrary time unit')
            plt.ylabel('Population numbers')
            string_for_title='SIROnGraphAgentBased, region num = '+str(n)+', k_SI ='+str(params['k_SI']) +', k_IR = '+ str(params['k_IR']) + ', S_ini = '+ str(params['S_ini']) + ', I_ini = '+str( params['I_ini']) + ', R_ini = '+ str(params['R_ini'])
            plt.title(string_for_title)
            string_for_savefig = 'SIROnGraphAgentBased_graph_type_'+str(params['graph_type'])+'_region_num='+str(n)+'_k_SI='+str(params['k_SI']) +'_k_IR='+ str(params['k_IR']) + '_S_ini='+ str(params['S_ini']) + '_I_ini='+str( params['I_ini']) + '_R_ini='+ str(params['R_ini'])
            plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
            plt.close()
                            
    if ('AllGraphs' in params['plot_options']):
        print("Running: plot option: 'AllGraphs'")
        ### --- loop over all saved graph and plot at each saved time step    
        for j in range(0,nbTimeSteps):
            print("saving ",j," th time step = " ,mytime[j])        
            # dic_x = outputs[j][1]
            dic_x = copy.deepcopy(outputs[j][1])
            dic_positions = nx.get_node_attributes(G,'positions')        
            dic_S = {}
            dic_I = {}
            dic_R = {}
            for i in dic_x.keys():
                dic_S[i] = len(dic_x[i][0])
                dic_I[i] = len(dic_x[i][1])
                dic_R[i] = len(dic_x[i][2])
            
            list_S_values = [v for v in dic_S.values()]
            list_I_values = [v for v in dic_I.values()]
            list_R_values = [v for v in dic_R.values()]                        
            fig20 = plt.figure(figsize=(16, 10)) # in inches!  
            string_for_title='Susceptible population on Graph at time = '+str(round(mytime[j],6))
            plt.title(string_for_title)
            ### --- drawing nodes and edges separately so we can capture collection for colobar
            ec = nx.draw_networkx_edges(G, pos=dic_positions, alpha=0.2)
            nc = nx.draw_networkx_nodes(G, pos=dic_positions, nodelist=G.nodes(), node_color=list_S_values, with_labels=False, node_size=100 ,cmap=plt.cm.viridis)#, cmap=plt.cm.viridis)
            plt.colorbar(nc)
            plt.axis('off')    
            string_for_savefig = 'Susceptible_population_on_Graph_graph_type='+str(params['graph_type'])+'_'+str('%08d'%j)
            plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')        
            plt.close()
            
            fig30 = plt.figure(figsize=(16, 10)) # in inches!
            string_for_title='Infected population on Graph at time = '+str(round(mytime[j],6))
            plt.title(string_for_title)
            ### --- drawing nodes and edges separately so we can capture collection for colobar
            ec = nx.draw_networkx_edges(G, pos=dic_positions, alpha=0.2)
            nc = nx.draw_networkx_nodes(G, pos=dic_positions, nodelist=G.nodes(), node_color=list_I_values, with_labels=False, node_size=100 ,cmap=plt.cm.viridis)#, cmap=plt.cm.viridis)
            plt.colorbar(nc)
            plt.axis('off')    
            string_for_savefig = 'Infected_population_on_Graph_graph_type='+str(params['graph_type'])+'_'+str('%08d'%j)
            plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
            plt.close()
            
            fig40 = plt.figure(figsize=(16, 10)) # in inches!            
            string_for_title='Recovered population on Graph at time = '+str(round(mytime[j],6))
            plt.title(string_for_title)
            ### --- drawing nodes and edges separately so we can capture collection for colobar
            ec = nx.draw_networkx_edges(G, pos=dic_positions, alpha=0.2)
            nc = nx.draw_networkx_nodes(G, pos=dic_positions, nodelist=G.nodes(), node_color=list_R_values, with_labels=False, node_size=100 ,cmap=plt.cm.viridis)#, cmap=plt.cm.viridis)
            plt.colorbar(nc)
            plt.axis('off')
            string_for_savefig = 'Recovered_population_on_Graph_graph_type='+str(params['graph_type'])+'_'+str('%08d'%j)
            plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
            plt.close()
    if ('FinalGraph' in params['plot_options']):
        print("Running: plot option: 'FinalGraph'")
        # ### --- plot graph with final species update == state in Graph
        dic_positions = nx.get_node_attributes(G,'positions')        
        dic_x = nx.get_node_attributes(G,'x')
        dic_S = {}
        dic_I = {}
        dic_R = {}
        for i in dic_x.keys():
            dic_S[i] = len(dic_x[i][0])
            dic_I[i] = len(dic_x[i][1])
            dic_R[i] = len(dic_x[i][2])
        
        list_S_values = [v for v in dic_S.values()]
        list_I_values = [v for v in dic_I.values()]
        list_R_values = [v for v in dic_R.values()]
        fig20 = plt.figure(figsize=(16, 10)) # in inches!  
        string_for_title='Susceptible population on Graph at time = '+str(round(G.graph['time'],6))
        plt.title(string_for_title)
        # ### --- drawing nodes and edges separately so we can capture collection for colobar
        ec = nx.draw_networkx_edges(G, pos=dic_positions, alpha=0.2)
        nc = nx.draw_networkx_nodes(G, pos=dic_positions, nodelist=G.nodes(), node_color=list_S_values, with_labels=False, node_size=100 ,cmap=plt.cm.viridis)#, cmap=plt.cm.viridis)
        plt.colorbar(nc)
        plt.axis('off')    
        string_for_savefig = 'Susceptible_population_on_Final_Graph_at_time='+str(round(G.graph['time'],6))
        plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
        plt.close()
        
        fig30 = plt.figure(figsize=(16, 10)) # in inches!
        string_for_title='Infected population on Graph at time = '+str(round(G.graph['time'],6))
        plt.title(string_for_title)
        # ### --- drawing nodes and edges separately so we can capture collection for colobar
        ec = nx.draw_networkx_edges(G, pos=dic_positions, alpha=0.2)
        nc = nx.draw_networkx_nodes(G, pos=dic_positions, nodelist=G.nodes(), node_color=list_I_values, with_labels=False, node_size=100 ,cmap=plt.cm.viridis)#, cmap=plt.cm.viridis)
        plt.colorbar(nc)
        plt.axis('off')    
        string_for_savefig = 'Infected_population_on_Final_Graph_at_time='+str(round(G.graph['time'],6))
        plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
        plt.close()
        
        fig40 = plt.figure(figsize=(16, 10)) # in inches!            
        string_for_title='Recovered population on Graph at time = '+str(round(G.graph['time'],6))
        plt.title(string_for_title)
        # ### --- drawing nodes and edges separately so we can capture collection for colobar
        ec = nx.draw_networkx_edges(G, pos=dic_positions, alpha=0.2)
        nc = nx.draw_networkx_nodes(G, pos=dic_positions, nodelist=G.nodes(), node_color=list_R_values, with_labels=False, node_size=100 ,cmap=plt.cm.viridis)#, cmap=plt.cm.viridis)
        plt.colorbar(nc)
        plt.axis('off')
        string_for_savefig = 'Recovered_population_on_Final_Graph_at_time='+str(round(G.graph['time'],6))
        plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')        
        plt.close()

    if ('TrackSuperSpreaders' in params['plot_options']):
        print("Running: plot option: 'TrackSuperSpreaders'")
        dic_x_ini = nx.get_node_attributes(G,'x_ini')
        nbNodes = G.number_of_nodes()
        nbStates = params['nbStates']
        nodes = G.nodes()
        all_nodes_all_x_over_time = {}        
        S_ini = params["S_ini"]
        I_ini = params["I_ini"]
        R_ini = params["R_ini"]
        N = S_ini + I_ini + R_ini
        # because of conservation of mass, there are the same numeber 
        # of individuals all along the simulation
        t0 = time.time()
        all_individuals_id = np.arange(0,N)
        
        dic_all_individuals_residence_nodes = {}
        for id_ in all_individuals_id:
            dic_all_individuals_residence_nodes[id_] = []
            for n in nodes:
                for j in range(0,nbStates):
                    if (id_ in dic_x_ini[n][j]):
                        starting_node = n
            ref_node = starting_node
            dic_all_individuals_residence_nodes[id_].append(ref_node)            
            for i in range(0,nbTimeSteps):
                dic_x = outputs[i][1]
                for n in nodes:
                    for j in range(0,nbStates):
                        if (id_ in dic_x[n][j]):
                            new_node = n
                if (new_node != ref_node):
                    # print('individual ', id_, ' initially in node ',  starting_node, ' hoped from node ', ref_node, ' to node ' , new_node, ' at time ' , mytime[i])
                    ref_node = new_node
                    dic_all_individuals_residence_nodes[id_].append(ref_node)
                
        # print(dic_all_individuals_residence_nodes)
        number_of_jumps = []
        for id_ in dic_all_individuals_residence_nodes.keys():
            number_of_jumps.append(len(dic_all_individuals_residence_nodes[id_]))
        number_of_jumps = np.array(number_of_jumps)
        print(number_of_jumps)
        t1 = time.time()
        print("number of jumps ran for: ", t1 - t0, " seconds")
        
        
        t0 = time.time()
        dic_total_individual_jumps = {}
        
        # create dictionnary of all individuals in every node (node is the key)
        dic_all_individuals_in_nodes_old = {}        
        dic_x = outputs[0][1]
        for n in nodes:
            all_individuals_in_node_old = [dic_x[n][i][j] for i in range(0,len(dic_x[n])) for j in range(0,len(dic_x[n][i]))]                
            dic_all_individuals_in_nodes_old[n] = all_individuals_in_node_old
        for i in range(1,nbTimeSteps): # start from second time step
            dic_x = outputs[i][1]
            individuals_that_jumped_at_this_time_step = []
            for n in nodes:
                all_individuals_in_node_old = dic_all_individuals_in_nodes_old[n]
                # concatenate lists of individuals
                all_individuals_in_node_new = [dic_x[n][i][j] for i in range(0,len(dic_x[n])) for j in range(0,len(dic_x[n][i]))]                
                # set is like a mathematical set
                # test if there was a jump                
                if (set(all_individuals_in_node_new) != set(all_individuals_in_node_old)):
                    if (len(set(all_individuals_in_node_new)) > len(set(all_individuals_in_node_old))):
                        individual_that_jumped = set(all_individuals_in_node_new) - set(all_individuals_in_node_old)
                    else:
                        individual_that_jumped = set(all_individuals_in_node_old) - set(all_individuals_in_node_new)
                    if (len(individual_that_jumped) != 1): 
                        print('two individuals jumped at the same time. ERROR')
                        print(individual_that_jumped)
                        exit(0)
                    else:
                        # transform set into integer                        
                        individual_that_jumped = [i for i in individual_that_jumped][0]                        
                        
                        # initialize in case it is the first time we register a jump
                        if not(individual_that_jumped in dic_total_individual_jumps):
                            dic_total_individual_jumps[individual_that_jumped] = 0
                        
                        # prevent double jumps at the same time step:
                        # e.g. we would count that id 1 leaves node 1 and that it arrives in node 2
                        if not(individual_that_jumped in individuals_that_jumped_at_this_time_step):
                            dic_total_individual_jumps[individual_that_jumped] += 1
                            individuals_that_jumped_at_this_time_step.append(individual_that_jumped)
                        # update new into old
                        all_individuals_in_node_old = all_individuals_in_node_new
                        dic_all_individuals_in_nodes_old[n] = all_individuals_in_node_old
                                                
        number_of_jumps_alternative = np.zeros(N)
        for id_ in dic_total_individual_jumps.keys():
            number_of_jumps_alternative[id_] = dic_total_individual_jumps[id_]
        print('\n', number_of_jumps_alternative)
        t1 = time.time()
        print("number of jumps alternative ran for: ", t1 - t0, " seconds")
        # --- plot distribution of number of jumps per id
        plt.hist(number_of_jumps, bins=20, density=True, label='number of jumps distribution')
        plt.legend()
        plt.xlabel('Number of jumps')
        plt.ylabel('PDF')
        string_for_title='SIROnGraphAgentBased, k_SI ='+str(params['k_SI']) +', k_IR = '+ str(params['k_IR']) + ', S_ini = '+ str(params['S_ini']) + ', I_ini = '+str( params['I_ini']) + ', R_ini = '+ str(params['R_ini'])
        plt.title(string_for_title)
        string_for_savefig = 'SIROnGraphAgentBased_number_of_jumps_distribution_graph_type_'+str(params['graph_type'])+'_k_SI='+str(params['k_SI']) +'_k_IR='+ str(params['k_IR']) + '_S_ini='+ str(params['S_ini']) + '_I_ini='+str( params['I_ini']) + '_R_ini='+ str(params['R_ini'])
        plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
        plt.close()

        # --- plot distribution of number of jumps per id
        plt.hist(number_of_jumps_alternative, bins=20, density=True, label='number of jumps distribution alternative')
        plt.legend()
        plt.xlabel('Number of jumps')
        plt.ylabel('PDF')
        string_for_title='SIROnGraphAgentBased, k_SI ='+str(params['k_SI']) +', k_IR = '+ str(params['k_IR']) + ', S_ini = '+ str(params['S_ini']) + ', I_ini = '+str( params['I_ini']) + ', R_ini = '+ str(params['R_ini'])
        plt.title(string_for_title)
        string_for_savefig = 'SIROnGraphAgentBased_number_of_jumps_distribution_alternatives_graph_type_'+str(params['graph_type'])+'_k_SI='+str(params['k_SI']) +'_k_IR='+ str(params['k_IR']) + '_S_ini='+ str(params['S_ini']) + '_I_ini='+str( params['I_ini']) + '_R_ini='+ str(params['R_ini'])
        plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
        plt.close()
