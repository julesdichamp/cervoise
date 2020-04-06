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


def plotData(params,G, days_to_plot):
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
    paramsToFit = {}
    # test params from fit 
    paramsToFit['k_SI' ] = 1.41350287619901e-10
    paramsToFit['k_IR' ] = 0.0025
    paramsToFit['k_IL' ] = 0.0005025625589838139
    paramsToFit['k_LD' ] = 0.025596063067312054
    paramsToFit['L_ini'] = 510.1404424811643
    paramsToFit['I_ini'] = 10**6
    paramsToFit['alpha'] = 0.1
    params["k_SI"]  = float(paramsToFit["k_SI"])
    params["k_IR"]  = float(paramsToFit["k_IR"])
    params["k_IL"]  = float(paramsToFit["k_IL"])
    params["k_LD"]  = float(paramsToFit["k_LD"])
    params["L_ini"] = int(paramsToFit["L_ini"])
    params["I_ini"] = int(paramsToFit["I_ini"])
    params["alpha"] = float(paramsToFit["alpha"])
    
    # compute W_c_t0 = int_t^t0 alpha * k_SI * S(s) * I(s) ds
    # tmp_func[0] = params['alpha'] * params["k_SI"] * x[0:i,0] * x[0:i,1]    
    # dates_before_lockdown = ['2020-03-01', '2020-03-02', '2020-03-03', '2020-03-04', '2020-03-05', '2020-03-06', '2020-03-07', '2020-03-08', '2020-03-09', '2020-03-10', '2020-03-11', '2020-03-12', '2020-03-13', '2020-03-14', '2020-03-15', '2020-03-16', '2020-03-17', '2020-03-18']
    # for _date in dates_before_lockdown:
        # I = params['alpha'] * params['total_covid_data']['Confirmed']
        # R = params['total_covid_data']['Recovered']
        # D = params['total_covid_data']['Deaths']
        # S = N - I - D - R
        # S = params['total_pop'] - 
    # W_c[i] = np.trapz(tmp_func, _time[0:i] ) # integral over time of I until _time[i]
    
    [params , G] = InitModels(params, G)
    outputs = SolveOneCompartment(params)
    nbTimes = len(outputs)
    W_c=np.zeros(nbTimes)    
    x = np.zeros([nbTimes,params['nbStates']])
    _time = np.zeros(nbTimes)
    print(nbTimes)
    for i in range(0,nbTimes):
        _time[i] = outputs[i][0]
        x[i,:] = outputs[i][1]
        tmp_func = params['alpha'] * params["k_SI"] * x[0:i+1,0] * x[0:i+1,1]
        W_c[i] = np.trapz(tmp_func, _time[0:i+1] ) # integral over time of I until _time[i]    
    
    W_c_t0 = params['alpha'] * params["k_SI"] * x[0,0] * x[0,1]
    W_c += W_c_t0
    
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
    
    plt.show()
