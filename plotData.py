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
import lmfit
sns.set_context("talk")
from models import InitModels, SolveOnGraph, SolveOneCompartment


def plotDataRegions(params,G, days_to_plot):
    """
    plot various data    
    """
    dic_data_covid = nx.get_node_attributes(G, name='data_covid')    
    
    fig1 = plt.figure(figsize=(16, 10)) # in inches!
    
    
    for n in G.nodes():        
        days=[]
        y = []
        for day in days_to_plot:
            days.append(day[5:])
            y.append(dic_data_covid[n][day]['Confirmed'])
        plt.plot(days,y, 'green', label='region_num='+str(n))
    
    plt.xlabel('Days')
    plt.ylabel('Population numbers')
    plt.legend()
    plt.savefig(params['outputs_dir']+'/Confirmed_all_regions_vs_time'+'.png')
    plt.close()
    
    fig2 = plt.figure(figsize=(16, 10)) # in inches!
    for n in G.nodes():        
        days=[]
        y = []
        for day in days_to_plot:
            days.append(day[5:])
            y.append(dic_data_covid[n][day]['Deaths'])
        plt.plot(days,y, 'green', label='region_num='+str(n))
    
    plt.xlabel('Days')
    plt.ylabel('Population numbers')
    plt.legend()
    plt.savefig(params['outputs_dir']+'/Deaths_all_regions_vs_time'+'.png')
    plt.close()
    
    fig3 = plt.figure(figsize=(16, 10)) # in inches!
    for n in G.nodes():        
        days=[]
        y = []
        for day in days_to_plot:
            days.append(day[5:])
            y.append(dic_data_covid[n][day]['Recovered'])
        plt.plot(days,y, 'green', label='region_num='+str(n))
    
    plt.xlabel('Days')
    plt.ylabel('Population numbers')
    plt.legend()
    plt.savefig(params['outputs_dir']+'/Recovered_all_regions_vs_time'+'.png')
    plt.close()



def plotDataFrance(params,G, days_to_plot):
    """
    plot covid data for entire france    
    """
        
    cases_confirmed = []
    cases_deaths = []
    cases_recovered = []
    days = []
    for day in days_to_plot:
        days.append(day[5:])
        cases_confirmed.append(params['total_covid_data'][day]['Confirmed'])
        cases_deaths.append(params['total_covid_data'][day]['Deaths'])
        cases_recovered.append(params['total_covid_data'][day]['Recovered'])    
    
    fig = plt.figure(figsize=(16, 10))
    plt.plot(cases_confirmed, label='Confirmed')
    plt.plot(cases_deaths, label='Deaths')
    plt.plot(cases_recovered, label='Recovered')    
    plt.xlabel('Days')
    plt.ylabel('Population numbers')
    plt.legend()
    plt.savefig(params['outputs_dir']+'/Covid_data_france_vs_time'+'.png')    
    plt.close()

    

def plotDataAndResults(params, G, days_to_plot):
    """
    plot W_c(t) = total number of confirmed at time t
    plot Deaths(t) = total number of deaths at time t
    """    
    total_Deaths = []
    total_Recovered = []
    total_Confirmed = []
    for _date in days_to_plot:
        total_Deaths.append(params['total_covid_data'][_date]['Deaths'])
        total_Recovered.append(params['total_covid_data'][_date]['Recovered'])
        total_Confirmed.append(params['total_covid_data'][_date]['Confirmed'])
    # plot data
    plt.figure('Deaths', figsize=(16, 10))
    p_deaths = plt.plot(days_to_plot, total_Deaths, linestyle=None, marker='o')
    plt.figure('Recovered', figsize=(16, 10))
    p_recovered = plt.plot(days_to_plot, total_Recovered, linestyle=None, marker='o')
    plt.figure('Confirmed', figsize=(16, 10))
    p_confirmed = plt.plot(days_to_plot, total_Confirmed, linestyle=None, marker='o')    
    
    [params , G] = InitModels(params, G)
    outputs = SolveOneCompartment(params)
    nbTimes = len(outputs)
    W_c=np.zeros(nbTimes)
    x = np.zeros([nbTimes,params['nbStates']])
    _time = np.zeros(nbTimes)    
    for i in range(0,nbTimes):
        _time[i] = outputs[i][0]
        x[i,:] = outputs[i][1]
        tmp_func = params['alpha'] * params["k_SI"] * x[0:i+1,0] * x[0:i+1,1]
        W_c[i] = np.trapz(tmp_func, _time[0:i+1] ) # integral over time of I until _time[i]  
        
    # W_c_t0 = params['alpha'] * params["k_SI"] * x[0,0] * x[0,1]
    # W_c += W_c_t0
    W_c  += params['total_covid_data'][params['time_start_day']]['Confirmed']
    
    plt.figure('Deaths', figsize=(16, 10))
    p_deaths = plt.plot(_time, x[:,4], linestyle=None, marker='o')
    plt.figure('Recovered', figsize=(16, 10))
    p_recovered = plt.plot(_time, params['alpha'] * x[:,2], linestyle=None, marker='o')
    plt.figure('Confirmed', figsize=(16, 10))
    p_confirmed = plt.plot(_time, W_c, linestyle=None, marker='o')
    

    # save
    plt.figure('Confirmed', figsize=(16, 10))
    string_for_savefig = 'Confirmed'
    plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
    plt.figure('Deaths', figsize=(16, 10))
    string_for_savefig = 'Deaths'
    plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')    
    plt.figure('Recovered', figsize=(16, 10))
    string_for_savefig = 'Recovered'
    plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
    
def plotDataAndSEIRLDResults(params, G, outputs, days_to_plot):
    """
    plot W_c(t) = total number of confirmed at time t
    plot Deaths(t) = total number of deaths at time t
    """    
    total_Deaths = []
    total_Recovered = []
    total_Confirmed = []
    for _date in days_to_plot:
        total_Deaths.append(params['total_covid_data'][_date]['Deaths'])
        total_Recovered.append(params['total_covid_data'][_date]['Recovered'])
        total_Confirmed.append(params['total_covid_data'][_date]['Confirmed'])
    
    # get total population for rescale of k_SEs
    dic_total_pop  = nx.get_node_attributes(G, name='total_pop')
    dic_simu_total_pop = copy.deepcopy(dic_total_pop)
    N_pop = 0
    for n in G.nodes():            
            N_pop += dic_simu_total_pop[n]
    
    nbTimes = len(outputs)
    W_c =np.zeros(nbTimes)
    x = np.zeros([nbTimes,params['nbStates']])
    _time = np.zeros(nbTimes)    
    for i in range(0,nbTimes):
        _time[i] = outputs[i][0]
        x[i,:] = outputs[i][1]        
        tmp_func = params['alpha'] * params["k_SE"] * x[0:i+1,0] * x[0:i+1,2] / N_pop
        W_c[i] = np.trapz(tmp_func, _time[0:i+1] ) # integral over time of I until _time[i]  
        
    print('final x = ' , x[-1,:])    
    fig2 = plt.figure(figsize=(16, 10)) # in inches!
    p_deaths = plt.plot(_time, x[:,5], linestyle=None, label='Deaths')    
    p_recovered = plt.plot(_time, x[:,3], linestyle=None, label='Recovered')    
    p_confirmed = plt.plot(_time, x[:,2], linestyle=None, label='Infected')    
    plt.legend()
    R0 = params['k_SE'] / (params['k_IR'] + params['k_IL'])
    string_for_savefig = 'RealSEIRLDWellMixed_TEST_k_SE='+str(params['k_SE'])+'_k_EI='+str(params['k_EI'])+'_k_IR='+str(params['k_IR'])+'_k_IL='+str(params['k_IL'])+'_k_LR='+str(params['k_LR'])+'_k_LD='+str(params['k_LD'])
    plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')    
    plt.show()
    
    W_c  += params['total_covid_data'][params['time_start_day']]['Confirmed']
    
    # plot data
    plt.figure('Deaths', figsize=(16, 10))
    p_deaths = plt.plot(days_to_plot, total_Deaths, linestyle=None, marker='o', label='data')
    plt.figure('Recovered', figsize=(16, 10))
    p_recovered = plt.plot(days_to_plot, total_Recovered, linestyle=None, marker='o', label='data')
    plt.figure('Confirmed', figsize=(16, 10))
    p_confirmed = plt.plot(days_to_plot, total_Confirmed, linestyle=None, marker='o', label='data')        
    
        
    print('final x = ' , x[-1,:])
    plt.figure('Deaths', figsize=(16, 10))
    p_deaths = plt.plot(_time, x[:,5], linestyle=None, label='model')
    plt.figure('Recovered', figsize=(16, 10))
    p_recovered = plt.plot(_time, params['alpha'] * x[:,3], linestyle=None, label='model')
    plt.figure('Confirmed', figsize=(16, 10))  
    p_confirmed = plt.plot(_time, W_c, linestyle=None, label='model')
    
    # save
    plt.figure('Confirmed', figsize=(16, 10))
    plt.title('France Confirmed cases model vs data')
    plt.legend()
    string_for_savefig = 'Confirmed'
    plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
    plt.figure('Deaths', figsize=(16, 10))
    plt.title('France Deaths cases model vs data')
    plt.legend()
    string_for_savefig = 'Deaths'
    plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')    
    plt.figure('Recovered', figsize=(16, 10))
    plt.title('France Recovered cases model vs data')
    plt.legend()
    string_for_savefig = 'Recovered'
    plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')    
    plt.show()


def plotDataAndSIRLDResults(params, G, outputs, days_to_plot):
    """
    plot W_c(t) = total number of confirmed at time t
    plot Deaths(t) = total number of deaths at time t
    """    
    total_Deaths = []
    total_Recovered = []
    total_Confirmed = []
    for _date in days_to_plot:
        total_Deaths.append(params['total_covid_data'][_date]['Deaths'])
        total_Recovered.append(params['total_covid_data'][_date]['Recovered'])
        total_Confirmed.append(params['total_covid_data'][_date]['Confirmed'])    
    
    # get total population for rescale of k_SEs
    dic_total_pop  = nx.get_node_attributes(G, name='total_pop')
    dic_simu_total_pop = copy.deepcopy(dic_total_pop)
    N_pop = 0
    for n in G.nodes():            
            N_pop += dic_simu_total_pop[n]
    
    nbTimes = len(outputs)
    W_c =np.zeros(nbTimes)
    x = np.zeros([nbTimes,params['nbStates']])
    _time = np.zeros(nbTimes)    
    for i in range(0,nbTimes):
        _time[i] = outputs[i][0]
        x[i,:] = outputs[i][1]        
        tmp_func = params['alpha'] * params["k_SI"] * x[0:i+1,0] * x[0:i+1,1] / N_pop
        W_c[i] = np.trapz(tmp_func, _time[0:i+1] ) # integral over time of I until _time[i]  
        
    print('final x = ' , x[-1,:])    
    fig2 = plt.figure(figsize=(16, 10)) # in inches!
    p_deaths = plt.plot(_time, x[:,4], linestyle=None, label='Deaths')    
    p_recovered = plt.plot(_time, x[:,2], linestyle=None, label='Recovered')    
    p_confirmed = plt.plot(_time, x[:,1], linestyle=None, label='Infected')    
    plt.legend()
    R0 = params['k_SI'] / (params['k_IR'] + params['k_IL'])
    string_for_savefig = 'RealSIRLDWellMixed_TEST_k_SE='+str(params['k_SI'])+'_k_IR='+str(params['k_IR'])+'_k_IL='+str(params['k_IL'])
    plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')    
    plt.show()


def plotDataAndSIRLDmodidifedResults(params, G, outputs, days_to_plot):
    """
    plot W_c(t) = total number of confirmed at time t
    plot Deaths(t) = total number of deaths at time t
    """    
    total_Deaths = []
    total_Recovered = []
    total_Confirmed = []
    for _date in days_to_plot:
        total_Deaths.append(params['total_covid_data'][_date]['Deaths'])
        total_Recovered.append(params['total_covid_data'][_date]['Recovered'])
        total_Confirmed.append(params['total_covid_data'][_date]['Confirmed'])    
    
    # get total population for rescale of k_SEs
    dic_total_pop  = nx.get_node_attributes(G, name='total_pop')
    dic_simu_total_pop = copy.deepcopy(dic_total_pop)
    N_pop = 0
    for n in G.nodes():            
            N_pop += dic_simu_total_pop[n]
    
    nbTimes = len(outputs)
    W_c =np.zeros(nbTimes)
    x = np.zeros([nbTimes,params['nbStates']])
    _time = np.zeros(nbTimes)    
    for i in range(0,nbTimes):
        _time[i] = outputs[i][0]
        x[i,:] = outputs[i][1]        
        tmp_func = params['alpha'] * params["k_SI"] * x[0:i+1,0] * x[0:i+1,1] / N_pop
        W_c[i] = np.trapz(tmp_func, _time[0:i+1] ) # integral over time of I until _time[i]  
        
    print('final x = ' , x[-1,:])    
    fig2 = plt.figure(figsize=(16, 10)) # in inches!
    p_deaths = plt.plot(_time, x[:,4], linestyle=None, label='Deaths')    
    p_recovered = plt.plot(_time, x[:,2], linestyle=None, label='Recovered')    
    p_confirmed = plt.plot(_time, x[:,1], linestyle=None, label='Infected')    
    plt.legend()
    R0 = params['k_SI'] / (params['k_IR'] + params['k_IL'])
    string_for_savefig = 'RealSIRLDmodifiedWellMixed_TEST_k_SI='+str(params['k_SI'])+'_k_IR='+str(params['k_IR'])+'_k_IL='+str(params['k_IL'])+'_k_LR='+str(params['k_LR'])
    plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')    
    plt.show()


def plotDataAndSEIRLDOnGraphResults(params, G, outputs, days_to_plot):
    """
    plot W_c(t) = total number of confirmed at time t
    plot Deaths(t) = total number of deaths at time t
    """       
    # save outputs
    f = open(params['outputs_dir']+'/outputs.txt',"w")
    f.write(str(outputs))
    f.close()
    
    # save params inputs
    f = open(params['outputs_dir']+'/params.txt',"w")
    f.write(str(params))
    f.close()
    
    G.graph['time'] = outputs[len(outputs)-1][0] # intialize time attribtue in graph
    nbTimes = len(outputs)
    print("nbTimes = " , nbTimes)
    nbStates = params['nbStates']
    dic_x = {}    
    time = np.zeros(nbTimes) 
    x = np.zeros([nbTimes, nbStates])    
    for i in range(0,nbTimes):        
        time[i] = outputs[i][0]
        dic_x = outputs[i][1]
        for n in dic_x.keys():
            for k in range(0, nbStates):
                x[i,k] += dic_x[n][k]      
    
    if ('WellMixed' in params['plot_options']):
        print("Running: plot option: 'WellMixed'")
        ### --- plot wellmixed model ( = sum over all compartments i.e. nodes of the graph)
        fig10 = plt.figure(figsize=(16, 10)) # in inches!
        # plt.plot(time,x[:,0], 'red', label='Susceptible')
        plt.plot(time,x[:,2], 'green', label='Infected')
        plt.plot(time,x[:,3], 'blue', label='Recovered')
        plt.plot(time,x[:,5], 'black', label='Deaths')
        plt.legend()
        plt.xlabel('Arbitrary time unit')
        plt.ylabel('Population numbers')
        string_for_title='SIROnGraph well mixed limit'#, k_SI ='+str(params['k_SI']) +', k_IR = '+ str(params['k_IR']) + ', S_ini = '+ str(params['S_ini']) + ', I_ini = '+str( params['I_ini']) + ', R_ini = '+ str(params['R_ini'])
        plt.title(string_for_title)
        string_for_savefig = 'SIROnGraph_WellMixed_limit'#_graph_type_'+str(params['graph_type'])+'_k_SI='+str(params['k_SI']) +'_k_IR='+ str(params['k_IR']) + '_S_ini='+ str(params['S_ini']) + '_I_ini='+str( params['I_ini']) + '_R_ini='+ str(params['R_ini'])
        plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
        plt.close()
        # plt.show()

    if ('AllNodes' in params['plot_options']):
        print("Running: plot option: 'AllNodes'")
        ### --- loop over all saved data and plot SIRWellMixed on each region
        # use dic_x available in node attributes from graph to get 
        # number of nodes and number of reactions                
        dic_names = nx.get_node_attributes(G,'names')
        nbReactions = params['nbReactions']
        nodes = G.nodes()
        nbNodes = len(nodes)
        all_nodes_all_x_over_time = {}
        for n in nodes:
            x = np.zeros([nbTimes, nbReactions])
            for i in range(0,nbTimes):
                dic_x = outputs[i][1]            
                for j in range(0,nbReactions):
                    x[i,j] = dic_x[n][j]                    
            all_nodes_all_x_over_time[n] = x
            
            plt.figure('AllRegionsInfected', figsize=(16, 10))
            p_all_regions_infected = plt.plot(time,x[:,2], label=str(dic_names[n]), linewidth = 3)
            
            fig10 = plt.figure('individual region', figsize=(16, 10)) # in inches!
            plt.plot(time,x[:,2], 'green', label='Infected')
            plt.plot(time,x[:,3], 'blue', label='Recovered')
            plt.plot(time,x[:,5], 'black', label='Deaths')
            plt.legend()
            plt.xlabel('Arbitrary time unit')
            plt.ylabel('Population numbers')
            string_for_title='SIROnGraph, region = '+str(dic_names[n])#+', k_SI ='+str(params['k_SI']) +', k_IR = '+ str(params['k_IR']) + ', S_ini = '+ str(params['S_ini']) + ', I_ini = '+str( params['I_ini']) + ', R_ini = '+ str(params['R_ini'])
            plt.title(string_for_title)
            string_for_savefig = 'SIROnGraph_region_num='+str(n)#+'_k_SI='+str(params['k_SI']) +'_k_IR='+ str(params['k_IR']) + '_S_ini='+ str(params['S_ini']) + '_I_ini='+str( params['I_ini']) + '_R_ini='+ str(params['R_ini'])
            plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
            plt.close()
        
        _fontsize = 30
        plt.figure('AllRegionsInfected', figsize=(16, 10))
        # plt.legend(fontsize=20, loc=1)
        plt.xlabel('Days', fontsize = _fontsize)
        plt.ylabel('Number of infected', fontsize = _fontsize)
        plt.xticks(fontsize=_fontsize)
        plt.yticks(fontsize=_fontsize)
        plt.tick_params(axis='both', which='major', length = 10, width = 2.5)
        plt.tick_params(axis='both', which='minor', length = 5, width = 1.5)
        plt.yscale('log')
        # plt.xlim([np.min(time), 1.8 * np.max(time)])
        string_for_savefig = 'SIROnGraph_all_regions_infected'
        plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
        plt.close()
        
    if 'FinalGraph' in params['plot_options']:
        print("Running: plot option: 'FinalGraph'")
        # ### --- plot graph with final species update == state in Graph        
        dic_positions = nx.get_node_attributes(G,'positions')        
        G.graph['time'] = outputs[len(outputs)-1][0] # intialize time attribtue in graph
        dic_x = outputs[len(outputs) - 1][1]
        nx.set_node_attributes(G, name='x', values=dic_x)
        nx.write_gpickle(G, params['outputs_dir']+'/Graph.gpickle', protocol = 2)
        dic_S = {}
        dic_I = {}
        dic_R = {}
        dic_D = {}
        for i in dic_x.keys():
            dic_S[i] = dic_x[i][0]
            dic_I[i] = dic_x[i][2]
            dic_R[i] = dic_x[i][3]
            dic_D[i] = dic_x[i][5]
        
        list_S_values = [v for v in dic_S.values()]
        list_I_values = [v for v in dic_I.values()]
        list_R_values = [v for v in dic_R.values()]
        list_D_values = [v for v in dic_D.values()]
        # fig20 = plt.figure(figsize=(16, 10)) # in inches!  
        # string_for_title='Susceptible population on Graph at time = '+str(round(G.graph['time'],6))
        # plt.title(string_for_title)
        ### --- drawing nodes and edges separately so we can capture collection for colobar
        # ec = nx.draw_networkx_edges(G, pos=dic_positions, alpha=0.2)
        # nc = nx.draw_networkx_nodes(G, pos=dic_positions, nodelist=G.nodes(), node_color=list_S_values, with_labels=False, node_size=100 ,cmap=plt.cm.viridis)#, cmap=plt.cm.viridis)
        # plt.colorbar(nc)
        # plt.axis('off')    
        # string_for_savefig = 'Susceptible_population_on_Final_Graph_at_time='+str(round(G.graph['time'],6))
        # plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
        # plt.close()
        
        fig30 = plt.figure(figsize=(16, 10)) # in inches!
        string_for_title='Infected population on Graph at time = '+str(round(G.graph['time'],6))
        plt.title(string_for_title)
        # ### --- drawing nodes and edges separately so we can capture collection for colobar
        ec = nx.draw_networkx_edges(G, pos=dic_positions, alpha=0.2)
        # nc = nx.draw_networkx_nodes(G, pos=dic_positions, nodelist=G.nodes(), node_color=list_I_values, with_labels=False, node_size=100 ,cmap=plt.cm.viridis)#, cmap=plt.cm.viridis)
        nc = nx.draw_networkx_nodes(G, pos=dic_positions, nodelist=G.nodes(), node_color=np.log(list_I_values), with_labels=False, node_size=100 ,cmap=plt.cm.viridis)#, cmap=plt.cm.viridis)
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
        # nc = nx.draw_networkx_nodes(G, pos=dic_positions, nodelist=G.nodes(), node_color=list_R_values, with_labels=False, node_size=100 ,cmap=plt.cm.viridis)#, cmap=plt.cm.viridis)
        nc = nx.draw_networkx_nodes(G, pos=dic_positions, nodelist=G.nodes(), node_color=np.log(list_R_values), with_labels=False, node_size=100 ,cmap=plt.cm.viridis)#, cmap=plt.cm.viridis)
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
        # nc = nx.draw_networkx_nodes(G, pos=dic_positions, nodelist=G.nodes(), node_color=list_D_values, with_labels=False, node_size=100 ,cmap=plt.cm.viridis)#, cmap=plt.cm.viridis)
        nc = nx.draw_networkx_nodes(G, pos=dic_positions, nodelist=G.nodes(), node_color=np.log(list_D_values), with_labels=False, node_size=100 ,cmap=plt.cm.viridis)#, cmap=plt.cm.viridis)
        plt.colorbar(nc)
        plt.axis('off')
        string_for_savefig = 'Deaths_population_on_Final_Graph_at_time='+str(round(G.graph['time'],6))
        plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')        
        plt.close()
