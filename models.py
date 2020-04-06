#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import geopandas as gpd
import pandas as pd
import numpy as np
import networkx as nx
import collections
import copy
from decimal import *
import fastrand

def randZeroOne(x):
    return fastrand.pcg32bounded(1000)/1000


def checkConservation(x_ini, x):
    
    if (isinstance(x_ini, collections.Mapping) and isinstance(x, collections.Mapping)):
        dic_x_ini = x_ini
        dic_x = x     
        total_mass_ini = 0.           
        total_mass = 0.
        for n in dic_x_ini.keys():
            for i in range(0,len(dic_x_ini[n])):
                total_mass_ini += dic_x_ini[n][i]
                total_mass += dic_x[n][i]
    else: # assume x_ini and x are numpy arrays a vector
        total_mass= np.sum(x)
        total_mass_ini= np.sum(x_ini)
    
    if abs(total_mass - total_mass_ini) > 0.0000001:
        print("conservation of mass is not valid")
        print("x_ini, total_mass_ini = " , x_ini, total_mass_ini)
        print("x, total_mass = " , x, total_mass)
        print("Exit here.")
        exit(-1)

def InitModels(params, G):
    """
    initialize model with options in 'params' dictionnary
    and add parameters to this dictionnary
    """    
    print('InitModels with: ', params["model_name"])
    """
    simple SIR model
    S + I -> 2I (with rate k_SI)
    I     -> R  (with rate k_IR)
    """     
    if (params["model_name"] == "SIRWellMixed"):
      
        k_SI = params["k_SI"]  # A = 1
        k_IR = params["k_IR"]  # A = 2
        S_ini = params["S_ini"]
        I_ini = params["I_ini"]
        R_ini = params["R_ini"]

        params['x_ini'] = np.array([S_ini, I_ini, R_ini])
        params['k_plus'] = np.array([k_SI, k_IR])
        # indices S = 0, I = 1, R = 2
        # N_ij = stoechiomatrix coefficients on left hand side of reactions
        # M_ij = stoechiomatrix coefficients on right hand side of reactions
        # i = state, j = reaction
        nbReactions = 2
        nbStates = 3
        N = np.zeros([nbStates, nbReactions])
        M = np.zeros([nbStates, nbReactions])
        
        #first reaction: S + I -> 2I
        N[0][0] = 1 ;  M[0][0] = 0
        N[1][0] = 1 ;  M[1][0] = 2
        N[2][0] = 0 ;  M[2][0] = 0
        #second reaction: I -> R
        N[0][1] = 0 ;  M[0][1] = 0
        N[1][1] = 1 ;  M[1][1] = 0
        N[2][1] = 0 ;  M[2][1] = 1
        
        R = M - N    
        params["nbReactions"] = nbReactions
        params["nbStates"] = nbStates
        params['N_matrix'] = N
        params['M_matrix'] = M
        params['R_matrix'] = R

    """
    simple SIRD model on well mixed hypothesis (D is death)
    S + I -> 2I (with rate k_SI)
    I     -> R  (with rate k_IR)
    I     -> D  (with rate k_ID)
    """
    if (params["model_name"] == "SIRDWellMixed"):
        k_SI = params["k_SI"]  # A = 1
        k_IR = params["k_IR"]  # A = 2
        k_ID = params["k_ID"]  # A = 3
        S_ini = params["S_ini"]
        I_ini = params["I_ini"]
        R_ini = params["R_ini"]
        D_ini = params["R_ini"]

        params['x_ini'] = np.array([S_ini, I_ini, R_ini, D_ini])
        params['k_plus'] = np.array([k_SI, k_IR, k_ID])
        # indices S = 0, I = 1, R = 2
        # N_ij = stoechiomatrix coefficients on left hand side of reactions
        # M_ij = stoechiomatrix coefficients on right hand side of reactions
        # i = state, j = reaction
        nbReactions = 3
        nbStates = 4
        N = np.zeros([nbStates, nbReactions])
        M = np.zeros([nbStates, nbReactions])
        
        # first reaction: S + I -> 2I
        N[0][0] = 1 ; M[0][0] = 0 
        N[1][0] = 1 ; M[1][0] = 2
        N[2][0] = 0 ; M[2][0] = 0
        N[3][0] = 0 ; M[3][0] = 0
        #second reaction: I -> R
        N[0][1] = 0 ; M[0][1] = 0
        N[1][1] = 1 ; M[1][1] = 0
        N[2][1] = 0 ; M[2][1] = 1
        N[3][1] = 0 ; M[3][1] = 0
        #third reaction: I -> D
        N[0][2] = 0 ; M[0][2] = 0
        N[1][2] = 1 ; M[1][2] = 0
        N[2][2] = 0 ; M[2][2] = 0
        N[3][2] = 0 ; M[3][2] = 1
        
        R = M - N
        params["nbReactions"] = nbReactions
        params["nbStates"] = nbStates
        params['N_matrix'] = N
        params['M_matrix'] = M
        params['R_matrix'] = R
        
    """
    SIR model with well mixed hypothesis and erlang waiting time
    for each reaction
    S + I   -> I + I_1 (with rate m * k_{SI})
    I_1     -> I_2     (with rate m * k_{SI})
    ...
    I_{m-1} -> I       (with rate m * k_{SI})
    
    I       -> R_1     (with rate m * k_{IR})
    ...
    R_{m-1} -> R       (with rate m * k_{IR})
    
    where m is the number of intermediate states    
    """ 
    if (params["model_name"] == "SIRerlangWellMixed"):
        k_SI = params["k_SI"]
        k_IR = params["k_IR"]
        S_ini = params["S_ini"]
        I_ini = params["I_ini"]
        R_ini = params["R_ini"]    
        m = params["nb_intermediate_states"]
        
        Lambda_SI = m * k_SI # recovery after 1 event
        Lambda_IR = m * k_IR # recovery after 1 event
            
        # indices S = 0, I_1 = 1, I_2 = 2... I (=I_m) = m
        # R_1 = m+1, ..., R (=R_m) = m+m
        # N_ij = stoechiomatrix coefficients on left hand side of reactions
        # M_ij = stoechiomatrix coefficients on right hand side of reactions
        # i = state, j = reaction
        nbReactions = 2 * m 
        # S, m * I_{i=1:m}, m*R_{i=1:m}
        # I_m = I, R_m = R
        nbStates = 1 + 2 * m
        
        # initialize population with known values
        x_ini = np.zeros(nbStates)
        x_ini[0] = S_ini
        x_ini[m] = I_ini
        x_ini[nbStates - 1] = R_ini        
        params['x_ini'] = x_ini
        
        # first m reactions  = Lambda_SI
        # second m reactions = Lambda_IR
        k_plus = np.zeros(nbReactions)
        k_plus[0:m] = Lambda_SI
        k_plus[m:nbReactions] = Lambda_IR    
        params['k_plus'] = k_plus
        
        N = np.zeros([nbStates, nbReactions])
        M = np.zeros([nbStates, nbReactions])
        
        # first reaction: S + I -> I + I1 (and all intermediate reactions)
        N[0][0] = 1 ; M[1][0] = 1
        N[m][0] = 1 ; M[m][0] += 1 # add 1 in case m = 1 and we need to have 2 on the right hand side
        # reactions I_1 -> I_2
        #           I_2 -> I_3
        #               ...
        #       I_{m-1} -> I_m = I
        # i ranges from 1 to m-1
        for i in range(1,m):
            # state indice is one more in right hand side than in left hand side
            N[i][i] = 1 #diagonal matrix
            M[i+1][i] = 1
        # reaction: I -> R_1
        N[m][m] = 1 ; M[m+1][m] = 1
        # reaections: R_1 -> R_2
        #                 ...
        #         R_{m-1} -> R_m = R
        # states number go from m + 1 to nbStates
        for i in range(m+1, nbReactions):
            N[i][i] = 1 # diagonal matrix
            M[i+1][i] = 1

        R = M - N  
        params["nbReactions"] = nbReactions
        params["nbStates"] = nbStates
        params['N_matrix'] = N
        params['M_matrix'] = M
        params['R_matrix'] = R
    
    """
    simple SIR model and hopping between nodes of graph        
    S + I -> 2I (with rate k_SI)
    I     -> R  (with rate k_IR)
    """        
    if (params["model_name"] == "SIROnGraph"):
        ### --- 1/ initializations
        # load initial state, and reactions values
        k_SI = params["k_SI"]  # A = 1
        k_IR = params["k_IR"]  # A = 2
        S_ini = params["S_ini"]
        I_ini = params["I_ini"]
        R_ini = params["R_ini"]

        x_ini = np.array([S_ini, I_ini, R_ini])
        k_plus = np.array([k_SI, k_IR])
        # indices S = 0, I = 1, R = 2
        # N_ij = stoechiomatrix coefficients on left hand side of reactions
        # M_ij = stoechiomatrix coefficients on right hand side of reactions
        # i = state, j = reaction
        nbReactions = 2
        nbStates = 3
        N = np.zeros([nbStates, nbReactions])
        M = np.zeros([nbStates, nbReactions])
        
        # first reaction: S + I -> 2I
        N[0][0] = 1 ;  M[0][0] = 0
        N[1][0] = 1 ;  M[1][0] = 2
        N[2][0] = 0 ;  M[2][0] = 0
        #second reaction: I -> R
        N[0][1] = 0 ;  M[0][1] = 0
        N[1][1] = 1 ;  M[1][1] = 0
        N[2][1] = 0 ;  M[2][1] = 1
                
        R = M - N
        
        params["nbReactions"] = nbReactions
        params["nbStates"] = nbStates
        params['N_matrix'] = N
        params['M_matrix'] = M
        params['R_matrix'] = R        
        
        dic_x_ini = {}
        dic_x = {}
        dic_kinetic = {}
        for n in G.nodes(): # loop over all compartments
            dic_x_ini[n] = x_ini/(G.number_of_nodes())
            dic_x[n] = x_ini/(G.number_of_nodes())
            dic_kinetic[n] = k_plus
        
        nx.set_node_attributes(G, name='x_ini', values=dic_x_ini)
        nx.set_node_attributes(G, name='x', values=dic_x)
        nx.set_node_attributes(G, name='k_plus', values=dic_kinetic)    
        
        # set edges attributes here
        dic_hopping_rates = {}
        for e in G.edges():
            dic_hopping_rates[e] = 0.01 # constant for now because TEST
        nx.set_edge_attributes(G, name='hopping_rates', values=dic_hopping_rates)
        
    """
    simple SIRD model with hopping between compartments      
    S + I -> 2I (with rate k_SI)
    I     -> R  (with rate k_IR)
    I     -> D  (with rate k_ID)
    """
    if (params["model_name"] == "SIRDOnGraph"):
        # load initial state, and reactions values
        k_SI = params["k_SI"]  # A = 1
        k_IR = params["k_IR"]  # A = 2
        k_ID = params["k_ID"]  # A = 2
        S_ini = params["S_ini"]        
        I_ini = params["I_ini"]
        R_ini = params["R_ini"]
        D_ini = params["D_ini"]
        
        x_ini = np.array([S_ini, I_ini, R_ini, D_ini])
        k_plus = np.array([k_SI, k_IR, k_ID])
        # indices S = 0, I = 1, R = 2
        # N_ij = stoechiomatrix coefficients on left hand side of reactions
        # M_ij = stoechiomatrix coefficients on right hand side of reactions
        # i = state, j = reaction
        nbReactions = 3
        nbStates = 4
        N = np.zeros([nbStates, nbReactions])
        M = np.zeros([nbStates, nbReactions])
        
        # first reaction: S + I -> 2I
        N[0][0] = 1 ; M[0][0] = 0 
        N[1][0] = 1 ; M[1][0] = 2
        N[2][0] = 0 ; M[2][0] = 0
        N[3][0] = 0 ; M[3][0] = 0
        #second reaction: I -> R
        N[0][1] = 0 ; M[0][1] = 0
        N[1][1] = 1 ; M[1][1] = 0
        N[2][1] = 0 ; M[2][1] = 1
        N[3][1] = 0 ; M[3][1] = 0
        #third reaction: I -> D
        N[0][2] = 0 ; M[0][2] = 0
        N[1][2] = 1 ; M[1][2] = 0
        N[2][2] = 0 ; M[2][2] = 0
        N[3][2] = 0 ; M[3][2] = 1
        
        R = M - N
        
        params["nbReactions"] = nbReactions
        params["nbStates"] = nbStates
        params['N_matrix'] = N
        params['M_matrix'] = M
        params['R_matrix'] = R    
        
        # JD: test here to set attributes to graph nodes
        dic_x_ini = {}
        dic_x = {}
        dic_kinetic = {}
        for n in G.nodes(): # loop over all compartments
            dic_x_ini[n] = x_ini/(G.number_of_nodes())
            dic_x[n] = x_ini/(G.number_of_nodes())
            dic_kinetic[n] = k_plus
        
        nx.set_node_attributes(G, name='x_ini', values=dic_x_ini)
        nx.set_node_attributes(G, name='x', values=dic_x)
        nx.set_node_attributes(G, name='k_plus', values=dic_kinetic)    
        
        # set edges attributes here
        dic_hopping_rates = {}
        for e in G.edges():
            dic_hopping_rates[e] = 0.01 # constant for now because TEST
        nx.set_edge_attributes(G, name='hopping_rates', values=dic_hopping_rates)    
    
    """
    SIR model with hopping between compartments and erlang waiting time
    for each reaction
    S + I   -> I + I_1 (with rate m * k_{SI})
    I_1     -> I_2     (with rate m * k_{SI})
    ...
    I_{m-1} -> I       (with rate m * k_{SI})
    
    I       -> R_1     (with rate m * k_{IR})
    ...
    R_{m-1} -> R       (with rate m * k_{IR})
    
    where m is the number of intermediate states    
    """    
    if (params["model_name"] == "SIRerlangOnGraph"):
        # load initial state, and reactions values
        k_SI = params["k_SI"]
        k_IR = params["k_IR"]
        S_ini = params["S_ini"]
        I_ini = params["I_ini"]
        R_ini = params["R_ini"]    
        m = params["nb_intermediate_states"]
        
        Lambda_SI = m * k_SI # recovery after 1 event
        Lambda_IR = m * k_IR # recovery after 1 event
            
        # indices S = 0, I_1 = 1, I_2 = 2... I (=I_m) = m
        # R_1 = m+1, ..., R (=R_m) = m+m
        # N_ij = stoechiomatrix coefficients on left hand side of reactions
        # M_ij = stoechiomatrix coefficients on right hand side of reactions
        # i = state, j = reaction
        nbReactions = 2 * m 
        # S, m * I_{i=1:m}, m*R_{i=1:m}
        # I_m = I, R_m = R
        nbStates = 1 + 2 * m
        
        # initialize population with known values
        x_ini = np.zeros(nbStates)
        x_ini[0] = S_ini
        x_ini[m] = I_ini
        x_ini[nbStates - 1] = R_ini
        
        # first m reactions  = Lambda_SI
        # second m reactions = Lambda_IR
        k_plus = np.zeros(nbReactions)
        k_plus[0:m] = Lambda_SI
        k_plus[m:nbReactions] = Lambda_IR    
        
        N = np.zeros([nbStates, nbReactions])
        M = np.zeros([nbStates, nbReactions])
        
        # first reaction: S + I -> I + I1 (and all intermediate reactions)
        N[0][0] = 1 ; M[1][0] = 1
        N[m][0] = 1 ; M[m][0] = 1
        # reactions I_1 -> I_2
        #           I_2 -> I_3
        #               ...
        #       I_{m-1} -> I_m = I
        # i ranges from 1 to m-1
        for i in range(1,m):
            # state indice is one more in right hand side than in left hand side
            N[i][i] = 1 #diagonal matrix
            M[i+1][i] = 1
        # reaction: I -> R_1
        N[m][m] = 1 ; M[m+1][m] = 1
        # reaections: R_1 -> R_2
        #                 ...
        #         R_{m-1} -> R_m = R
        # states number go from m + 1 to nbStates
        for i in range(m+1, nbReactions):
            N[i][i] = 1 # diagonal matrix
            M[i+1][i] = 1
            
        R = M - N    
        
        params["nbReactions"] = nbReactions
        params["nbStates"] = nbStates
        params['N_matrix'] = N
        params['M_matrix'] = M
        params['R_matrix'] = R   
        
        # JD: test here to set attributes to graph nodes
        dic_x_ini = {}
        dic_x = {}
        dic_kinetic = {}
        for n in G.nodes(): # loop over all compartments
            dic_x_ini[n] = x_ini/(G.number_of_nodes())
            dic_x[n] = x_ini/(G.number_of_nodes())
            dic_kinetic[n] = k_plus
        
        nx.set_node_attributes(G, name='x_ini', values=dic_x_ini)
        nx.set_node_attributes(G, name='x', values=dic_x)
        nx.set_node_attributes(G, name='k_plus', values=dic_kinetic)    
        
        # set edges attributes here
        dic_hopping_rates = {}
        for e in G.edges():
            dic_hopping_rates[e] = 0.01 # constant for now because TEST
        nx.set_edge_attributes(G, name='hopping_rates', values=dic_hopping_rates)
    
    """
    temporary name for testing stuff 
    """
    if (params["model_name"] == "RealSIRDOnGraph"):
        # load initial state, and reactions values
        R_0 = params["R_0"]  # A = 1
        k_IR = params["k_IR"]  # A = 2
        k_ID = params["k_ID"]  # A = 2
                
        # indices S = 0, I = 1, R = 2
        # N_ij = stoechiomatrix coefficients on left hand side of reactions
        # M_ij = stoechiomatrix coefficients on right hand side of reactions
        # i = state, j = reaction
        nbReactions = 3
        nbStates = 4
        N = np.zeros([nbStates, nbReactions])
        M = np.zeros([nbStates, nbReactions])
        
        # first reaction: S + I -> 2I
        N[0][0] = 1 ; M[0][0] = 0 
        N[1][0] = 1 ; M[1][0] = 2
        N[2][0] = 0 ; M[2][0] = 0
        N[3][0] = 0 ; M[3][0] = 0
        #second reaction: I -> R
        N[0][1] = 0 ; M[0][1] = 0
        N[1][1] = 1 ; M[1][1] = 0
        N[2][1] = 0 ; M[2][1] = 1
        N[3][1] = 0 ; M[3][1] = 0
        #third reaction: I -> D
        N[0][2] = 0 ; M[0][2] = 0
        N[1][2] = 1 ; M[1][2] = 0
        N[2][2] = 0 ; M[2][2] = 0
        N[3][2] = 0 ; M[3][2] = 1
        
        R = M - N
        
        params["nbReactions"] = nbReactions
        params["nbStates"] = nbStates
        params['N_matrix'] = N
        params['M_matrix'] = M
        params['R_matrix'] = R    
        
        dic_data_covid = nx.get_node_attributes(G, name='data_covid')
        dic_total_pop  = nx.get_node_attributes(G, name='total_pop')
        dic_simu_total_pop = copy.deepcopy(dic_total_pop)
        dic_simu_data_covid = copy.deepcopy(dic_data_covid)        
        dic_kinetic = {}
        dic_x_ini = {}
        dic_x = {}
        for n in G.nodes():
            # take 1% of total population and round to floor to make sure 
            # we work with integer population numbers            
            dic_simu_total_pop[n] = dic_simu_total_pop[n]
            # get fraction of confirmed cases             
            I = dic_simu_data_covid[n][params["time_start_day"]]['Confirmed']
            # get fraction of deaths
            D = dic_simu_data_covid[n][params["time_start_day"]]['Deaths']
            R = dic_simu_data_covid[n][params["time_start_day"]]['Recovered']
            S = dic_simu_total_pop[n] - I - D - R
            
            # adjust k_SI to each node population
            k_SI = (k_IR + k_ID) * R_0 / dic_simu_total_pop[n]            
            k_plus = np.array([k_SI, k_IR, k_ID])
            # add all kinetic to dictionnary
            dic_kinetic[n] = k_plus
            dic_x_ini[n] = np.array([S,I,R,D])
            dic_x[n] = np.array([S,I,R,D])
            print("node, x_ini = " , n , dic_x_ini[n])

        nx.set_node_attributes(G, name='x_ini', values=dic_x_ini)
        nx.set_node_attributes(G, name='x', values=dic_x)
        nx.set_node_attributes(G, name='k_plus', values=dic_kinetic)
        nx.set_node_attributes(G, name='simu_total_pop', values=dic_simu_total_pop)
        
        # set edges attributes here
        dic_hopping_rates = {}
        for e in G.edges():
            dic_hopping_rates[e] = 0.  # do not consider hopping after lockdown
        nx.set_edge_attributes(G, name='hopping_rates', values=dic_hopping_rates) 
    
    """
    temporary name for testing stuff 
    """
    if (params["model_name"] == "RealSIRDWellMixed"):
        # load initial state, and reactions values
        R_0 = params["R_0"]  # A = 1
        k_IR = params["k_IR"]  # A = 2
        k_ID = params["k_ID"]  # A = 2
                
        # indices S = 0, I = 1, R = 2
        # N_ij = stoechiomatrix coefficients on left hand side of reactions
        # M_ij = stoechiomatrix coefficients on right hand side of reactions
        # i = state, j = reaction
        nbReactions = 3
        nbStates = 4
        N = np.zeros([nbStates, nbReactions])
        M = np.zeros([nbStates, nbReactions])
        
        # first reaction: S + I -> 2I
        N[0][0] = 1 ; M[0][0] = 0 
        N[1][0] = 1 ; M[1][0] = 2
        N[2][0] = 0 ; M[2][0] = 0
        N[3][0] = 0 ; M[3][0] = 0
        #second reaction: I -> R
        N[0][1] = 0 ; M[0][1] = 0
        N[1][1] = 1 ; M[1][1] = 0
        N[2][1] = 0 ; M[2][1] = 1
        N[3][1] = 0 ; M[3][1] = 0
        #third reaction: I -> D
        N[0][2] = 0 ; M[0][2] = 0
        N[1][2] = 1 ; M[1][2] = 0
        N[2][2] = 0 ; M[2][2] = 0
        N[3][2] = 0 ; M[3][2] = 1
        
        R = M - N
        
        params["nbReactions"] = nbReactions
        params["nbStates"] = nbStates
        params['N_matrix'] = N
        params['M_matrix'] = M
        params['R_matrix'] = R    
        
        dic_data_covid = nx.get_node_attributes(G, name='data_covid')
        dic_total_pop  = nx.get_node_attributes(G, name='total_pop')
        dic_simu_total_pop = copy.deepcopy(dic_total_pop)
        dic_simu_data_covid = copy.deepcopy(dic_data_covid)  
        S = 0 ; I = 0 ; R = 0 ; D = 0;
        # get data for entire france
        for n in G.nodes():
            # take 1% of total population and round to floor to make sure 
            # we work with integer population numbers            
            dic_simu_total_pop[n] = dic_simu_total_pop[n]
            # get fraction of confirmed cases             
            I += dic_simu_data_covid[n][params["time_start_day"]]['Confirmed']
            # get fraction of deaths
            D += dic_simu_data_covid[n][params["time_start_day"]]['Deaths']
            R += dic_simu_data_covid[n][params["time_start_day"]]['Recovered']
            S += dic_simu_total_pop[n] - I - D - R

        # get N = total population
        N = S + I + R + D

        # adjust k_SI to each node population
        k_SI = (k_IR + k_ID) * R_0 / N
        k_plus = np.array([k_SI, k_IR, k_ID])
        x_ini = np.array([S, I, R, D])
        params['k_plus'] = k_plus
        params['x_ini'] = x_ini
    
    """
    model from marc lavielle 
    http://webpopix.org/covidix19.html
    In our framework it becomes:
    S + I --> 2I
    I     --> R
    I     --> L
    L     --> D
    """
    if (params["model_name"] == "RealSIRLDWellMixed"):
        # load initial state, and reactions values
        k_SI = params["k_SI"]  # A = 0
        k_IR = params["k_IR"]  # A = 1
        k_IL = params["k_IL"]  # A = 2
        k_LD = params["k_LD"]  # A = 3
        # indices S = 0, I = 1, R = 2
        # N_ij = stoechiomatrix coefficients on left hand side of reactions
        # M_ij = stoechiomatrix coefficients on right hand side of reactions
        # i = state, j = reaction
        nbReactions = 4
        nbStates = 5
        N = np.zeros([nbStates, nbReactions])
        M = np.zeros([nbStates, nbReactions])
        
        # first reaction: S + I -> 2I
        N[0][0] = 1 ; M[0][0] = 0 
        N[1][0] = 1 ; M[1][0] = 2
        N[2][0] = 0 ; M[2][0] = 0
        N[3][0] = 0 ; M[3][0] = 0
        N[4][0] = 0 ; M[4][0] = 0
        #second reaction: I -> R
        N[0][1] = 0 ; M[0][1] = 0
        N[1][1] = 1 ; M[1][1] = 0
        N[2][1] = 0 ; M[2][1] = 1
        N[3][1] = 0 ; M[3][1] = 0
        N[4][1] = 0 ; M[4][1] = 0
        #third reaction: I -> L
        N[0][2] = 0 ; M[0][2] = 0
        N[1][2] = 1 ; M[1][2] = 0
        N[2][2] = 0 ; M[2][2] = 0
        N[3][2] = 0 ; M[3][2] = 1
        N[4][2] = 0 ; M[4][2] = 0
        #fourth reaction: L -> D
        N[0][3] = 0 ; M[0][3] = 0
        N[1][3] = 0 ; M[1][3] = 0
        N[2][3] = 0 ; M[2][3] = 0
        N[3][3] = 1 ; M[3][3] = 0
        N[4][3] = 0 ; M[4][3] = 1
        
        R = M - N
        
        params["nbReactions"] = nbReactions
        params["nbStates"] = nbStates
        params['N_matrix'] = N
        params['M_matrix'] = M
        params['R_matrix'] = R    
        
        dic_data_covid = nx.get_node_attributes(G, name='data_covid')
        dic_total_pop  = nx.get_node_attributes(G, name='total_pop')
        dic_simu_total_pop = copy.deepcopy(dic_total_pop)
        dic_simu_data_covid = copy.deepcopy(dic_data_covid)  
        S = 0 ; I = 0 ; R = 0 ; D = 0; N = 0;
        # get data for entire france
        for n in G.nodes():                        
            # I_confirmed = alpha * I, thus I = 1/alpha * I_confirmed
            # I += dic_simu_data_covid[n][params["time_start_day"]]['Confirmed']
            # get fraction of deaths
            D += dic_simu_data_covid[n][params["time_start_day"]]['Deaths']
            R += dic_simu_data_covid[n][params["time_start_day"]]['Recovered']            
            dic_simu_total_pop[n]

            N += dic_simu_total_pop[n]
        
        # get N = total population        
        L = params['L_ini'] # people who are isolated and will die 
        # I = np.floor((1/params['alpha']) * I)
        I = params['I_ini']
        S = N - I - D - R
        
        # adjust k_SI to each node population        
        k_plus = np.array([k_SI, k_IR, k_IL, k_LD])
        x_ini = np.array([S, I, R, L, D])
        params['k_plus'] = k_plus
        params['x_ini'] = x_ini
        
        
    return params, G
        
def SolveOneCompartment(params):
    """
    solve master equation with one compartment    
    """
    nbReactions = params["nbReactions"]
    nbStates = params["nbStates"]
    k_plus         = params['k_plus']
    N              = params['N_matrix']
    R              = params['R_matrix']
    x_ini          = params['x_ini']
    
    x = x_ini
    t = params["time_start"]
    outputs = []  # we will fill this list with data to be plotted
    count_loops = 0
    while t < params["time_end"]:        
        outputs.append((t, x))
        if ((count_loops % 100000) == 0):
            print("t, count_loops, x = " , t, count_loops, x)
        # forward reaction rates
        t_forward_plus = np.zeros(nbReactions)
        for j in range(0, nbReactions):
            t_forward_plus[j] = 1.0
            for i in range(0, nbStates):
                # y is the product from x[i] * (x[i] - 1) * ... * (x[i] - N[i,j] + 1)
                y = 1.0
                # only considers states that are present in current reaction
                # see Example 1 from DD powerpoint
                if (abs(N[i,j]) > 0.00000001):
                    for k in range(0,int(N[i,j])):                        
                        y = y * (x[i] -k)
                        
                t_forward_plus[j] = t_forward_plus[j] * y
            
            # finally multiply the forward reaction rate by the reaction rate k+
            t_forward_plus[j] = t_forward_plus[j] * k_plus[j]
                        
        # total rate which the state can be left with
        W = np.sum(t_forward_plus)
        xi = randZeroOne()
        # if the total rate is 0, nothing more can happen anymore
        if (abs(W) <= 0.00000001):
            print("W = " , W)
            break
        
        # Gillepsie time of next event
        dt = -1 / W * np.log(1 - xi)
        t = t + dt

        # sum reactions forward rates to choose reaction
        next_states = np.zeros(nbReactions)
        next_states[0] = t_forward_plus[0]
        for j in range(1, nbReactions):
            next_states[j] = next_states[j - 1] + t_forward_plus[j]
        next_states = next_states / W

        # choose reaction
        xi = randZeroOne()
        if xi <= next_states[0]:
            id_reaction = 0
        else:
            for j in range(1, nbReactions):
                if (next_states[j - 1] < xi) and (xi < next_states[j]):
                    id_reaction = j

        # set new state according to reaction        
        x = x + R[:, id_reaction] 
        
        count_loops += 1
        
    return outputs

def SolveOnGraph(params, G):
    """
    solve master equation with several compartments (nodes of graph)
    with hopping processes
    """
    nbReactions = params["nbReactions"]
    nbStates = params["nbStates"]
    N              = params['N_matrix']
    R              = params['R_matrix']
    
    # intialize states and time    
    t = params["time_start"]
    G.graph['time'] = t # intialize time attribtue in graph
    outputs = []  # we will fill this list with data to be plotted
    count_loops = 0 # loop counter
    ### --- 2/ loop in time
    try: # so we can exit the loop at any time with ctrl + C
        while (t < params["time_end"]):
            dic_x = nx.get_node_attributes(G,'x')
            dic_k_plus = nx.get_node_attributes(G,'k_plus')
            
            # save time as an attribute in graph at each time step
            G.graph['time'] = t
                    
            # save every n time steps
            if ((count_loops % params['save_frequency']) == 0):            
                dic_x_for_save = copy.deepcopy(dic_x)
                x_print = np.zeros(nbStates)
                for key in dic_x_for_save.keys():
                    for ii in range(0,nbStates):
                        x_print[ii] += dic_x_for_save[key][ii]
                # print("t, count_loops, x = " , t, count_loops, x_print)
                if ((count_loops % 10000) == 0):
                    print("t, count_loops, x = " , t, count_loops, x_print)
                outputs.append([t, dic_x_for_save])
            
            # compute W_inside: loop on every node 
            W_inside = 0.0
            dic_t_forward_plus = {}        
            ### --- 2-1/ compute total inside rate = sum of all forward kinetic rates)
            for n in G.nodes(): # loop over all compartments
                x = dic_x[n]
                k_plus = dic_k_plus[n]
                # forward reaction rates
                t_forward_plus = np.zeros(nbReactions)
                for j in range(0, nbReactions):
                    t_forward_plus[j] = 1.0
                    for i in range(0, nbStates):
                        # y is the product from x[i] * (x[i] - 1) * ... * (x[i] - N[i,j] + 1)
                        y = 1.0
                        # only considers states that are present in current reaction
                        # (see Example 1 from DD powerpoint)
                        if (abs(N[i,j]) > 0.00000001):
                            for k in range(0,int(N[i,j])):
                                y = y * (x[i] -k)                                                
                        t_forward_plus[j] = t_forward_plus[j] * y
                    
                    # finally multiply the forward reaction rate by the reaction rate k+
                    t_forward_plus[j] = t_forward_plus[j] * k_plus[j]
                  
                dic_t_forward_plus[n] = t_forward_plus # save forward reaction rates for all nodes
                
                # total rate which the state can be left with
                W_node = np.sum(t_forward_plus)
                xi = randZeroOne()
                W_inside += W_node
            
            ### --- 2-2/ compute total between rate = total rate of leaving the state with a 
            # hopping process = sum of all hopping rates d_ij * xi
            dic_t_between = {}
            W_between = 0.0
            dic_hopping_rates = nx.get_edge_attributes(G,'hopping_rates')
            for ni in G.nodes(): # loop over all nodes
                for nj in G.neighbors(ni): # loop over all neighbors
                    # not sure if it can happen that a node sees itself as a neighbor but just in case...
                    if (ni != nj):
                        eij=(ni,nj) # get edge between the two nodes                                        
                        # complicated stuff to turn the edge the other way around
                        # in case it is not found in edge list
                        found_edge=False
                        for e in G.edges():
                            if (e == eij):
                                found_edge=True
                        if not(found_edge):
                            eij = (nj,ni)
                                            
                        tmp = np.zeros(nbStates)
                        for i in range(0, nbStates): # loop over all states                        
                            W_between += dic_x[ni][i] * dic_hopping_rates[eij]
                            tmp[i] = dic_x[ni][i] * dic_hopping_rates[eij]
                        #VERY IMPORTANT CHECK
                        dic_t_between[(ni,nj)] = tmp # save all between rates 
                    
            ### --- 2-3/ compute total rate and Gillepsie time step
            W_total = W_inside + W_between
            
            # Gillepsie time of next event using total rate 
            dt = -1 / W_total * np.log(1 - xi)
            t = t + dt            
            
            ### --- 2-4/ fill a vector that contains all sum of states: kinetics and hopping
            ###          there is probably a simpler way to do it but for now we play "by the book"
            ###          It allows to choose randomly but according to the weight of each reaction, 
            ###          the next event to trigger (hopping process of ONE state or ONE reaction)
            next_states_sum = np.zeros(nbReactions*len(G.nodes()) + nbStates*len(G.nodes()) * (len(G.nodes()) - 1))
            # fill states nbReactions*len(G.nodes())
            kk=0        
            for n in G.nodes(): # loop over all compartments            
                for j in range(0,nbReactions):
                    # initialization. We do it here because it is annoying to work
                    # with graph nodes indices
                    # it is only done for the first compartment and the first reaction
                    if (kk==0):
                        next_states_sum[kk] = dic_t_forward_plus[n][j]
                    else:
                        next_states_sum[kk] = next_states_sum[kk-1] + dic_t_forward_plus[n][j]
                    kk+=1
            # now fill states: nbStates*len(G.nodes()) * (len(G.nodes()) - 1)
            for ni in G.nodes(): # loop over all nodes
                for nj in G.neighbors(ni): # loop over all neighbors
                    if (ni != nj):                        
                        for i in range(0, nbStates): # loop over all states
                            next_states_sum[kk] = next_states_sum[kk-1] + dic_t_between[(ni,nj)][i]
                            kk+=1
                            
            # normalize so that all values are bteween 0 and 1
            next_states_sum = next_states_sum / W_total
            
            ### --- 2-5/ choose reaction according to a random uniform variable drop        
            xi = randZeroOne()
            if xi <= next_states_sum[0]:
                id_reaction = 0
            else:
                for j in range(1, len(next_states_sum)):
                    if (next_states_sum[j - 1] < xi) and (xi < next_states_sum[j]):
                        id_reaction = j

            ### --- 2-6/ Set new state
            ###          look for the id of the reaction and trigger the corresponding event
            ###          there is probably a simpler way to do it.
            if (id_reaction < nbReactions*len(G.nodes())): # we set a kinetic process
                # print("kinetic")
                # IMPORTANT: this assumes that nodes range from 0 to nb_nodes - 1
                node = np.floor(id_reaction / nbReactions)
                react = id_reaction%nbReactions
                # dic_x[node] = dic_x[node] - R[:, react]
                # below is the version in notes, we use this for consistency
                # also R has to be changed as R = M - N
                dic_x[node] = dic_x[node] + R[:, react]
                # #update attributes of nodes right away
                nx.set_node_attributes(G, name='x', values=dic_x)    
                if ((dic_x[node][0] < 0) or (dic_x[node][1] < 0) or (dic_x[node][2] < 0) ):
                    print('negative kinetic: id_reaction, node, react, dic_x[node] = ' , id_reaction, node, react, dic_x[node])
                    exit(-1)
            else: # we set a hopping process 
                # print("hopping")            
                kk = nbReactions*len(G.nodes())
                kk_init = kk
                for ni in G.nodes(): # loop over all nodes
                    for nj in G.neighbors(ni): # loop over all neighbors
                        if (ni != nj):
                            for i in range(0, nbStates): # loop over all states
                                if (kk==id_reaction):
                                    dic_x[ni][i] -= 1 # add one                            
                                    dic_x[nj][i] += 1 # remove one
                                    # update attributes right away
                                    nx.set_node_attributes(G, name='x', values=dic_x)
                                    if ((dic_x[ni][i] < 0) or (dic_x[nj][i] < 0)):
                                        print('negative hopping: id_reaction, ni, nj, i, kk, dic_x[ni][i], dic_x[nj][i]  = ' , id_reaction,  ni, nj, i, kk, dic_x[ni][i], dic_x[nj][i])
                                        print('dic_t_between[(',ni,',',nj,')][',i,' = ' , dic_t_between[(ni,nj)][i])
                                        print('dic_t_between[(',nj,',',ni,')][',i,' = ' , dic_t_between[(nj,ni)][i])
                                        print('dic_t_between = ' , dic_t_between)
                                        exit(-1)
                                kk+=1      
            ### check there is conservation
            ### -- comment it for now?
            # dic_x_ini = nx.get_node_attributes(G,'x_ini')
            # dic_x = nx.get_node_attributes(G,'x')
            # checkConservation(dic_x_ini, dic_x)
            count_loops+=1        
    except KeyboardInterrupt:
        pass
    return outputs , G

