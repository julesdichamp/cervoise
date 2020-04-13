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
    string_for_title='SIRerlangWellMixed, k_SI ='+str(params['k_SI'])+' , nb_intermediate_states = ' +str(params['nb_intermediate_states']) + ', S_ini = '+ str(params['S_ini']) + ', I_ini = '+str( params['I_ini'])
    plt.title(string_for_title)
    string_for_savefig = 'SIRerlangWellMixed_k_SI='+str(params['k_SI']) +'_nb_intermediate_states=' + str(params['nb_intermediate_states']) + '_S_ini='+ str(params['S_ini']) + '_I_ini='+str( params['I_ini'])
    plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
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
    string_for_title='SIRerlangOtherWellMixed, k_SI ='+str(params['k_SI'])+' , nb_intermediate_states = ' +str(params['nb_intermediate_states']) + ', S_ini = '+ str(params['S_ini']) + ', I_ini = '+str( params['I_ini'])
    plt.title(string_for_title)
    string_for_savefig = 'SIRerlangOtherWellMixed_k_SI='+str(params['k_SI']) +'_nb_intermediate_states=' + str(params['nb_intermediate_states']) + '_S_ini='+ str(params['S_ini']) + '_I_ini='+str( params['I_ini'])
    plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
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

