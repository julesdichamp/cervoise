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
    # plt.show()
    
    # fig_confirmed = plt.figure('confirmed',figsize=(16, 10))
    # plt.plot(cases_confirmed, label='Confirmed')    
    # plt.xlabel('Days')
    # plt.ylabel('Population numbers')
    # plt.legend()
    # plt.savefig(params['outputs_dir']+'/Covid_confirmed_france_vs_time'+'.png')    
    # plt.close()
    
    # fig_confirmed = plt.figure('deaths',figsize=(16, 10))
    # plt.plot(cases_deaths, label='Deaths')
    # plt.xlabel('Days')
    # plt.ylabel('Population numbers')
    # plt.legend()
    # plt.savefig(params['outputs_dir']+'/Covid_deaths_france_vs_time'+'.png')    
    # plt.close()
    
    # fig_confirmed = plt.figure('recovered',figsize=(16, 10))
    # plt.plot(cases_recovered, label='Recovered')    
    # plt.xlabel('Days')
    # plt.ylabel('Population numbers')
    # plt.legend()
    # plt.savefig(params['outputs_dir']+'/Covid_recovered_france_vs_time'+'.png')        
    # plt.close()
    

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
    
    # plt.show()

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
    # plot data
    plt.figure('Deaths', figsize=(16, 10))
    p_deaths = plt.plot(days_to_plot, total_Deaths, linestyle=None, marker='o')
    plt.figure('Recovered', figsize=(16, 10))
    p_recovered = plt.plot(days_to_plot, total_Recovered, linestyle=None, marker='o')
    plt.figure('Confirmed', figsize=(16, 10))
    p_confirmed = plt.plot(days_to_plot, total_Confirmed, linestyle=None, marker='o')    
        
    nbTimes = len(outputs)
    W_c =np.zeros(nbTimes)
    x = np.zeros([nbTimes,params['nbStates']])
    _time = np.zeros(nbTimes)    
    for i in range(0,nbTimes):
        _time[i] = outputs[i][0]
        x[i,:] = outputs[i][1]        
        tmp_func = params['alpha'] * params["k_SE"] * x[0:i+1,0] * x[0:i+1,2]
        W_c[i] = np.trapz(tmp_func, _time[0:i+1] ) # integral over time of I until _time[i]  
        
    W_c  += params['total_covid_data'][params['time_start_day']]['Confirmed']
        
    plt.figure('Deaths', figsize=(16, 10))
    p_deaths = plt.plot(_time, x[:,5], linestyle=None)
    plt.figure('Recovered', figsize=(16, 10))
    p_recovered = plt.plot(_time, params['alpha'] * x[:,3], linestyle=None)
    plt.figure('Confirmed', figsize=(16, 10))
    p_confirmed = plt.plot(_time, W_c, linestyle=None)
    
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
