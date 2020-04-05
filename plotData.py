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
    
    
    # fig2 = plt.figure(figsize=(16, 10)) # in inches!
    # for day in days_to_plot:
        # for n in G.nodes():
            # dic_data_covid[n][day]['Deaths']
    # fig3 = plt.figure(figsize=(16, 10)) # in inches!
    # for day in days_to_plot:
        # for n in G.nodes():
            # dic_data_covid[n][day]['Recovered']
            
    # plt.plot(time,x[:,1], 'green', label='Infected')
    # plt.plot(time,x[:,2], 'blue', label='Recovered')
    # plt.legend()
    # plt.xlabel('Arbitrary time unit')
    # plt.ylabel('Population numbers')
    # string_for_title='SIRWellMixed, k_SI ='+str(params['k_SI']) +', k_IR = '+ str(params['k_IR']) + ', S_ini = '+ str(params['S_ini']) + ', I_ini = '+str( params['I_ini']) + ', R_ini = '+ str(params['R_ini'])
    # plt.title(string_for_title)
    # string_for_savefig = 'SIRWellMixed_k_SI='+str(params['k_SI']) +'_k_IR='+ str(params['k_IR']) + '_S_ini='+ str(params['S_ini']) + '_I_ini='+str( params['I_ini']) + '_R_ini='+ str(params['R_ini'])
    # plt.savefig(params['outputs_dir']+'/'+str(string_for_savefig)+'.png')
