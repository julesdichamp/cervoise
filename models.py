#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import geopandas as gpd
import pandas as pd
import numpy as np
import networkx as nx
import collections
import copy

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
    else: # assume x_ini and x are numpy arrays a vector
        total_mass= np.sum(x)
        total_mass_ini= np.sum(x_ini)
    
    if abs(total_mass - total_mass_ini) > 0.0000001:
        print("conservation of mass is not valid")
        print("x_ini, total_mass_ini = " , x_ini, total_mass_ini)
        print("x, total_mass = " , x, total_mass)
        print("Exit here.")
        exit(-1)    
        
    
    

def SIRWellMixed(params):
    """
    simple SIR model on well mixed hypothesis:
    Susceptible -> Infectious
    Infectious -> Recovered
    k_SI is rate from Susceptible to Infectious
    k_IR is rate from Infectious to Recovered
    """
    print("SIRWellMixed")
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
    N[0][0] = 1
    N[1][1] = 1
    M[1][0] = 1
    M[2][1] = 1
    R = N - M


    # JD: this loop is quite general and can be used whatever the model
    # only minor modifiactions have to be made to use on another model
    
    # intialize states and time
    x = x_ini
    t = params["time_start"]
    outputs = []  # we will fill this list with data to be plotted
    while t < params["time_end"]:        
        print(t, x)
        outputs.append((t, x))
        # forward reaction rates
        t_forward_plus = np.zeros(nbReactions)
        for j in range(0, nbReactions):
            t_forward_plus[j] = 1.0
            for i in range(0, nbStates):
                # JD: very ineffective way to compute t_forward_plus: if x[i] 
                # is too big factorial will explode
                # t_forward_plus[j] = (
                    # t_forward_plus[j]
                    # * np.math.factorial(x[i])
                    # / np.math.factorial(x[i] - N[i, j])
                # )
                # y is the product from x[i] * (x[i] - 1) * ... * (x[i] - N[i,j] + 1)
                y = 1.0
                # only considers states that are present in current reaction
                # see Example 1 from DD powerpoint
                if (abs(N[i,j]) > 0.00000001):
                    for k in range(0,int(N[i,j])):                        
                        y = y * (x[i] -k)
                        
                t_forward_plus[j] = t_forward_plus[j] * y
            
            # finally multiply the forward reaction rate by the reaction rate k+
            t_forward_plus[j] = t_forward_plus[j] * k_plus[j]
                        
        # total rate which the state can be left with
        W = np.sum(t_forward_plus)
        xi = np.random.uniform(0, 1)
        # if the total rate is 0, nothing more can happen anymore
        if (abs(W) <= 0.00000001):
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
        xi = np.random.uniform(0, 1)        
        if xi <= next_states[0]:
            id_reaction = 0
        else:
            for j in range(1, nbReactions):
                if (next_states[j - 1] < xi) and (xi < next_states[j]):
                    id_reaction = j

        # set new state according to reaction
        x = x - R[:, id_reaction]

    return outputs

def SIROnGraph(params, G):
    """
    Compute SIR model on graph
    """
    print("SIROnGraph")
    print("WARNING: this only works if number of population per state is an integer")
    print("TODO: ")
    ### --- 1/ initializations
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
    N[0][0] = 1
    N[1][1] = 1
    M[1][0] = 1
    M[2][1] = 1
    R = N - M
    
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
    
    # set edges attributes here
    dic_hopping_rates = {}
    for e in G.edges():
        dic_hopping_rates[e] = 1.0 # constant for now because TEST
    nx.set_edge_attributes(G, name='hopping_rates', values=dic_hopping_rates)
    
    # intialize states and time    
    t = params["time_start"]
    G.graph['time'] = t # intialize time attribtue in graph
    outputs = []  # we will fill this list with data to be plotted
    count_loops = 0 # loop counter
    ### --- 2/ loop in time
    while t < params["time_end"]:            
        dic_x = nx.get_node_attributes(G,'x')
        dic_k_plus = nx.get_node_attributes(G,'k_plus')
        
        # save time as an attribute in graph at each time step
        G.graph['time'] = t
                
        # save every n time steps
        if ((count_loops % params['save_frequency']) == 0):
            print("t, count_loops = " , t, count_loops)
            dic_x_for_save = copy.deepcopy(dic_x)
            outputs.append([t, dic_x_for_save])
        
        # compute W_inside: loop on every node 
        W_inside = 0.0
        dic_t_forward_plus = {}        
        ### --- 2-1/ compute total inside rate = sum of all forward kinetic rates)
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
                    # only considers states that are present in current reaction
                    # (see Example 1 from DD powerpoint)
                    if (abs(N[i,j]) > 0.00000001):
                        for k in range(0,int(N[i,j])):
                            y = y * (x[i] -k)                                                
                    t_forward_plus[j] = t_forward_plus[j] * y
                
                # finally multiply the forward reaction rate by the reaction rate k+
                t_forward_plus[j] = t_forward_plus[j] * k_plus[j]
              
            dic_t_forward_plus[n] = t_forward_plus # save forward reaction rates for all nodes
            
            # total rate which the state can be left with
            W_node = np.sum(t_forward_plus)
            xi = np.random.uniform(0, 1)
            W_inside += W_node
        
        ### --- 2-2/ compute total between rate = total rate of leaving the state with a 
        # hopping process = sum of all hopping rates d_ij * xi
        dic_t_between = {}
        W_between = 0.0
        dic_hopping_rates = nx.get_edge_attributes(G,'hopping_rates')
        for ni in G.nodes(): # loop over all nodes
            for nj in G.neighbors(ni): # loop over all neighbors
                # not sure if it can happen that a node sees itself as a neighbor but just in case...
                if (ni != nj):
                    eij=(ni,nj) # get edge between the two nodes                                        
                    # complicated stuff to turn the edge the other way around
                    # in case it is not found in edge list
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
                    dic_t_between[(ni,nj)] = tmp # save all between rates 
                
        ### --- 2-3/ compute total rate and Gillepsie time step
        W_total = W_inside + W_between
        
        # Gillepsie time of next event using total rate 
        dt = -1 / W_total * np.log(1 - xi)
        t = t + dt        
        
        
        ### --- 2-4/ fill a vector that contains all sum of states: kinetics and hopping
        ###          there is probably a simpler way to do it but for now we play "by the book"
        ###          It allows to choose randomly but according to the weight of each reaction, 
        ###          the next event to trigger (hopping process of ONE state or ONE reaction)
        next_states_sum = np.zeros(nbReactions*len(G.nodes()) + nbStates*len(G.nodes()) * (len(G.nodes()) - 1))
        # fill states nbReactions*len(G.nodes())
        kk=0        
        for n in G.nodes(): # loop over all compartments            
            for j in range(0,nbReactions):
                # initialization. We do it here because it is annoying to work
                # with graph nodes indices
                # it is only done for the first compartment and the first reaction
                if (kk==0):
                    next_states_sum[kk] = dic_t_forward_plus[n][j]
                else:
                    next_states_sum[kk] = next_states_sum[kk-1] + dic_t_forward_plus[n][j]
                kk+=1
        # now fill states: nbStates*len(G.nodes()) * (len(G.nodes()) - 1)
        for ni in G.nodes(): # loop over all nodes
            for nj in G.neighbors(ni): # loop over all neighbors
                if (ni != nj):                        
                    for i in range(0, nbStates): # loop over all states
                        next_states_sum[kk] = next_states_sum[kk-1] + dic_t_between[(ni,nj)][i]
                        kk+=1
                        
        # normalize so that all values are bteween 0 and 1
        next_states_sum = next_states_sum / W_total
        
        ### --- 2-5/ choose reaction according to a random uniform variable drop        
        xi = np.random.uniform(0, 1)
        if xi <= next_states_sum[0]:
            id_reaction = 0
        else:
            for j in range(1, len(next_states_sum)):
                if (next_states_sum[j - 1] < xi) and (xi < next_states_sum[j]):
                    id_reaction = j
    
        ### --- 2-6/ Set new state
        ###          look for the id of the reaction and trigger the corresponding event
        ###          there is probably a simpler way to do it.
        if (id_reaction < nbReactions*len(G.nodes())): # we set a kinetic process
            # print("kinetic")
            # IMPORTANT: this assumes that nodes range from 0 to nb_nodes - 1
            node = np.floor(id_reaction / nbReactions)
            react = id_reaction%nbReactions
            dic_x[node] = dic_x[node] - R[:, react]
            # #update attributes of nodes right away
            nx.set_node_attributes(G, name='x', values=dic_x)    
            if ((dic_x[node][0] < 0) or (dic_x[node][1] < 0) or (dic_x[node][2] < 0) ):
                print('negative kinetic: id_reaction, node, react, dic_x[node] = ' , id_reaction, node, react, dic_x[node])
                exit(-1)
        else: # we set a hopping process 
            # print("hopping")            
            kk = nbReactions*len(G.nodes())
            kk_init = kk
            for ni in G.nodes(): # loop over all nodes
                for nj in G.neighbors(ni): # loop over all neighbors
                    if (ni != nj):
                        for i in range(0, nbStates): # loop over all states
                            if (kk==id_reaction):
                                dic_x[ni][i] -= 1 # add one                            
                                dic_x[nj][i] += 1 # remove one
                                # update attributes right away
                                nx.set_node_attributes(G, name='x', values=dic_x)
                                if ((dic_x[ni][i] < 0) or (dic_x[nj][i] < 0)):
                                    print('negative hopping: id_reaction, ni, nj, i, kk, dic_x[ni][i], dic_x[nj][i]  = ' , id_reaction,  ni, nj, i, kk, dic_x[ni][i], dic_x[nj][i])
                                    print('dic_t_between[(',ni,',',nj,')][',i,' = ' , dic_t_between[(ni,nj)][i])
                                    print('dic_t_between[(',nj,',',ni,')][',i,' = ' , dic_t_between[(nj,ni)][i])
                                    print('dic_t_between = ' , dic_t_between)
                                    exit(-1)
                            kk+=1      
        ### check there is conservation
        ### -- comment it for now?
        # dic_x_ini = nx.get_node_attributes(G,'x_ini')
        # dic_x = nx.get_node_attributes(G,'x')
        # checkConservation(dic_x_ini, dic_x)
        count_loops+=1        
            
    return outputs , G
    
def SIROnRegions(params, graph):
    """
    TODO
    """


def SIROnDepartments(params, graph):
    """
    TODO
    """
