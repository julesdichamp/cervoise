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

def randZeroOne():
    return fastrand.pcg32bounded(1000000)/1000000

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

        Lambda_SI = k_SI * m # recovery after 1 event
        Lambda_IR = k_IR * m # recovery after 1 event

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
        print("k_plus = ", k_plus)

    """
    SIR model with well mixed hypothesis and erlang waiting time
    for the first reaction reaction
    S + I   -> I + I_1 (with rate m * k_{SI})
    I_1     -> I_2     (with rate m * k_{SI})
    ...
    I_{m-1} -> I       (with rate m * k_{SI})

    I       -> R       (with rate k_{IR})    
    where m is the number of intermediate states
    """
    if (params["model_name"] == "SIRerlang1WellMixed"):
        
        k_SI = params["k_SI"]
        k_IR = params["k_IR"]
        S_ini = params["S_ini"]
        I_ini = params["I_ini"]
        R_ini = params["R_ini"]
        m = params["nb_intermediate_states"]

        Lambda_SI = k_SI * m # recovery after 1 event        

        # indices S = 0, I_1 = 1, I_2 = 2... I (=I_m) = m
        # R = indice m +1
        # N_ij = stoechiomatrix coefficients on left hand side of reactions
        # M_ij = stoechiomatrix coefficients on right hand side of reactions
        # i = state, j = reaction
        nbReactions = 1 + m
        # S, m * I_{i=1:m}, 1*R        
        nbStates = m + 2

        # initialize population with known values
        x_ini = np.zeros(nbStates)
        x_ini[0] = S_ini
        x_ini[m] = I_ini
        x_ini[nbStates - 1] = R_ini
        params['x_ini'] = x_ini

        # first m reactions  = Lambda_SI
        # second m reactions = Lambda_IR
        k_plus = np.zeros(nbReactions)
        # k_plus[0:m] = Lambda_SI
        # k_plus[m:nbReactions] = k_IR
        k_plus[0] = k_SI * (m)
        k_plus[1:m] = k_SI * (m) * S_ini
        k_plus[m:nbReactions] = k_IR
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
        # reaction: I -> R
        N[m][m] = 1 ; M[m+1][m] = 1

        R = M - N
        params["nbReactions"] = nbReactions
        params["nbStates"] = nbStates
        params['N_matrix'] = N
        params['M_matrix'] = M
        params['R_matrix'] = R
        print('N = \n' , N)
        print('M = \n' , M)
        print('R = \n' , R)
        print("k_plus = ", k_plus)    

    """
    SIR model with well mixed hypothesis and erlang waiting time
    for the second reaction reaction
    S + I   -> 2 I (with rate k_{SI})
    I -> R_1              (with rate m * k_{IR})
    R_1 -> R_2            (with rate m * k_{IR})
    ...
    R_{m-1} -> R  (= R_m) (with rate m * k_{IR})
    where m is the number of intermediate states
    """
    if (params["model_name"] == "SIRerlang2WellMixed"):
        
        k_SI = params["k_SI"]
        k_IR = params["k_IR"]
        S_ini = params["S_ini"]
        I_ini = params["I_ini"]
        R_ini = params["R_ini"]
        m = params["nb_intermediate_states"]

        Lambda_IR = k_IR * m # recovery after 1 event        
        
        # indices S = 0, I_1 = 1, I_2 = 2... I (=I_m) = m
        # R = indice m +1
        # N_ij = stoechiomatrix coefficients on left hand side of reactions
        # M_ij = stoechiomatrix coefficients on right hand side of reactions
        # i = state, j = reaction
        nbReactions = 1 + m
        # S, m * I_{i=1:m}, 1*R        
        nbStates = m + 2

        # initialize population with known values
        x_ini = np.zeros(nbStates)
        x_ini[0] = S_ini
        x_ini[1] = I_ini
        x_ini[nbStates - 1] = R_ini
        params['x_ini'] = x_ini

        # first m reactions  = Lambda_SI
        # second m reactions = Lambda_IR
        k_plus = np.zeros(nbReactions)
        k_plus[0] = k_SI
        k_plus[1:nbReactions] = k_IR * m
        params['k_plus'] = k_plus

        N = np.zeros([nbStates, nbReactions])
        M = np.zeros([nbStates, nbReactions])

        # first reaction: S + I -> I + I1 (and all intermediate reactions)
        N[0][0] = 1 ; M[1][0] = 1
        N[1][0] = 1 ; M[1][0] += 1 # add 1 in case m = 1 and we need to have 2 on the right hand side
        # reaction: I -> R_1
        N[1][1] = 1 ; M[1+1][1] = 1
        # reaections: R_1 -> R_2
        #                 ...
        #         R_{m-1} -> R_m = R
        # states number go from m + 1 to nbStates
        for i in range(1+1, nbReactions):
            N[i][i] = 1 # diagonal matrix
            M[i+1][i] = 1

        R = M - N
        params["nbReactions"] = nbReactions
        params["nbStates"] = nbStates
        params['N_matrix'] = N
        params['M_matrix'] = M
        params['R_matrix'] = R
        print('N = \n' , N)
        print('M = \n' , M)
        print('R = \n' , R)
        print("k_plus = ", k_plus)    

    """
    SI model with well mixed hypothesis a
    S + I   -> I + I (with rate k_{SI})    
    where m is the number of intermediate states
    """
    if (params["model_name"] == "SIWellMixed"):
        
        k_SI = params["k_SI"]        
        S_ini = params["S_ini"]
        I_ini = params["I_ini"]        

        # indices S = 0, I_1 = 1, I_2 = 2... I (=I_m) = m
        # R = indice m +1
        # N_ij = stoechiomatrix coefficients on left hand side of reactions
        # M_ij = stoechiomatrix coefficients on right hand side of reactions
        # i = state, j = reaction
        nbReactions = 1
        # S, m * I_{i=1:m}, 1*R        
        nbStates = 2

        # initialize population with known values
        x_ini = np.zeros(nbStates)
        x_ini[0] = S_ini
        x_ini[1] = I_ini        
        params['x_ini'] = x_ini

        # first m reactions  = Lambda_SI
        # second m reactions = Lambda_IR
        k_plus = np.zeros(nbReactions)
        k_plus[0] = k_SI        
        params['k_plus'] = k_plus

        N = np.zeros([nbStates, nbReactions])
        M = np.zeros([nbStates, nbReactions])

        # first reaction: S + I -> I + I1 (and all intermediate reactions)
        N[0][0] = 1 ; M[1][0] = 1
        N[1][0] = 1 ; M[1][0] += 1 # add 1 in case m = 1 and we need to have 2 on the right hand side

        R = M - N
        params["nbReactions"] = nbReactions
        params["nbStates"] = nbStates
        params['N_matrix'] = N
        params['M_matrix'] = M
        params['R_matrix'] = R
        print('N = \n' , N)
        print('M = \n' , M)
        print('R = \n' , R)
        print("k_plus = ", k_plus)    

    """
    SI model with well mixed hypothesis and erlang waiting time
    for the first reaction reaction
    S + I   -> I + I_1 (with rate m * k_{SI})
    I_1     -> I_2     (with rate m * k_{SI})
    ...
    I_{m-1} -> I       (with rate m * k_{SI})    
    where m is the number of intermediate states
    """
    if (params["model_name"] == "SIerlangWellMixed"):
        
        k_SI = params["k_SI"]
        S_ini = params["S_ini"]
        I_ini = params["I_ini"]
        m = params["nb_intermediate_states"]

        Lambda_SI = k_SI * m # recovery after 1 event        

        # indices S = 0, I_1 = 1, I_2 = 2... I (=I_m) = m
        # R = indice m +1
        # N_ij = stoechiomatrix coefficients on left hand side of reactions
        # M_ij = stoechiomatrix coefficients on right hand side of reactions
        # i = state, j = reaction
        nbReactions =  m
        # S, m * I_{i=1:m}, 1*R        
        nbStates = m + 1

        # initialize population with known values
        x_ini = np.zeros(nbStates)
        x_ini[0] = S_ini
        x_ini[m] = I_ini
        # x_ini[1] = I_ini
        params['x_ini'] = x_ini

        # first m reactions  = Lambda_SI
        # second m reactions = Lambda_IR
        k_plus = np.zeros(nbReactions)
        k_plus[0] = k_SI
        k_plus[1:m] = k_SI * m  * S_ini
        # k_plus[1:m] = k_SI * (m)
        # k_plus[0] = k_SI   * (m)
        # k_plus[1:m] = k_SI * (m)
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

        R = M - N
        params["nbReactions"] = nbReactions
        params["nbStates"] = nbStates
        params['N_matrix'] = N
        params['M_matrix'] = M
        params['R_matrix'] = R
        print('N = \n' , N)
        print('M = \n' , M)
        print('R = \n' , R)
        print("k_plus = ", k_plus) 
        print("x_ini = ", x_ini) 

    """
    SI model with well mixed hypothesis and erlang waiting time
    for the first reaction reaction and we modify it so that 
    that there is no I_m trigger. [USEFUL TO COMPUTE WAITING TIMES]
    S + I   -> I + I_1 (with rate m * k_{SI})
    I_1     -> I_2     (with rate m * k_{SI})
    ...
    I_{m-1} -> I_{m-1}       (with rate m * k_{SI})    
    where m is the number of intermediate states
    """
    if (params["model_name"] == "SIerlangWellMixedModified"):
        
        k_SI = params["k_SI"]
        S_ini = params["S_ini"]
        I_ini = params["I_ini"]
        m = params["nb_intermediate_states"]

        Lambda_SI = k_SI * m # recovery after 1 event        

        # indices S = 0, I_1 = 1, I_2 = 2... I (=I_m) = m
        # R = indice m +1
        # N_ij = stoechiomatrix coefficients on left hand side of reactions
        # M_ij = stoechiomatrix coefficients on right hand side of reactions
        # i = state, j = reaction
        nbReactions =  m
        # S, m * I_{i=1:m}, 1*R        
        nbStates = m + 1

        # initialize population with known values
        x_ini = np.zeros(nbStates)
        x_ini[0] = S_ini
        x_ini[m] = I_ini
        params['x_ini'] = x_ini
        
        k_plus = np.zeros(nbReactions)
        k_plus[0] =   k_SI * (m)
        k_plus[1:m] = k_SI * (m) #* S_ini
        # k_plus[0] = k_SI   * (m)
        # k_plus[1:m] = k_SI * (m)
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

        R = M - N
        params["nbReactions"] = nbReactions
        params["nbStates"] = nbStates
        params['N_matrix'] = N
        params['M_matrix'] = M
        params['R_matrix'] = R
        print('N = \n' , N)
        print('M = \n' , M)
        print('R = \n' , R)
        print("k_plus = ", k_plus) 
        print("x_ini = ", x_ini) 

    """
    SI model with erlang and fixed Susceptible population
    """
    if (params["model_name"] == "SIerlangFixSWellMixed"):
        
        k_SI = params["k_SI"]
        S_ini = params["S_ini"]
        I_ini = params["I_ini"]
        m = params["nb_intermediate_states"]

        Lambda_SI = k_SI * m # recovery after 1 event        

        # indices S = 0, I_1 = 1, I_2 = 2... I (=I_m) = m
        # R = indice m +1
        # N_ij = stoechiomatrix coefficients on left hand side of reactions
        # M_ij = stoechiomatrix coefficients on right hand side of reactions
        # i = state, j = reaction
        nbReactions =  m
        # S, m * I_{i=1:m}, 1*R        
        nbStates = m + 1

        # initialize population with known values
        x_ini = np.zeros(nbStates)
        x_ini[0] = S_ini
        x_ini[m] = I_ini
        params['x_ini'] = x_ini

        # first m reactions  = Lambda_SI
        # second m reactions = Lambda_IR
        k_plus = np.zeros(nbReactions)
        # k_plus[0] = k_SI
        # k_plus[1:m] = k_SI * S_ini
        k_plus[0] = k_SI * m
        k_plus[1:m] = k_SI * (m) * S_ini
        params['k_plus'] = k_plus

        N = np.zeros([nbStates, nbReactions])
        M = np.zeros([nbStates, nbReactions])

        # first reaction: S + I -> S + I + I1 (and all intermediate reactions)        
        # (S is an enzyme)
        N[0][0] = 1
        N[m][0] = 1
        M[0][0] = 1 # so that R[0][0] = 0 
        M[1][0] = 1
        M[m][0] += 1 # add 1 in case m = 1 and we need to have 2 on the right hand side
        # reactions I_1 -> I_2
        #           I_2 -> I_3
        #               ...
        #       I_{m-1} -> I_m = I
        # i ranges from 1 to m-1
        for i in range(1,m):
            # state indice is one more in right hand side than in left hand side
            N[i][i] = 1 #diagonal matrix
            M[i+1][i] = 1        

        R = M - N
        
        # small hack to ensure that Susceptible population is kept constant:        
        # R[0][0] = 0 # it should -1 here
        
        params["nbReactions"] = nbReactions
        params["nbStates"] = nbStates
        params['N_matrix'] = N
        params['M_matrix'] = M
        params['R_matrix'] = R
        print('N = \n' , N)
        print('M = \n' , M)
        print('R = \n' , R)
        print("k_plus = ", k_plus)
        print("x_ini = ", x_ini)

    """
    SI model with well mixed hypothesis and another type of
    erlang waiting time
    for the first reaction reaction    
    where m is the number of intermediate states
    """
    if (params["model_name"] == "SIerlangOtherWellMixed"):
        
        k_SI = params["k_SI"]
        S_ini = params["S_ini"]
        I_ini = params["I_ini"]
        m = params["nb_intermediate_states"]

        Lambda_SI = k_SI * m # recovery after 1 event        

        # indices S = 0, I_1 = 1, I_2 = 2... I (=I_m) = m
        # R = indice m +1
        # N_ij = stoechiomatrix coefficients on left hand side of reactions
        # M_ij = stoechiomatrix coefficients on right hand side of reactions
        # i = state, j = reaction
        nbReactions =  m + m-1
        # S, m * I_{i=1:m}, 1*R        
        nbStates = m + 1

        # initialize population with known values
        x_ini = np.zeros(nbStates)
        x_ini[0] = S_ini
        x_ini[1] = I_ini
        params['x_ini'] = x_ini

        # first m reactions  = Lambda_SI
        # second m reactions = Lambda_IR
        k_plus = np.zeros(nbReactions)
        k_plus[0:m] = k_SI * m
        k_plus[m:nbReactions] = k_SI * m * S_ini
        # k_plus[0:m] = k_SI * m
        # k_plus[m:nbReactions] = k_SI * m
        params['k_plus'] = k_plus

        N = np.zeros([nbStates, nbReactions])
        M = np.zeros([nbStates, nbReactions])
        
        # reactions S + I_[k} -> I_1 + I_k
        # there are m of them
        for i in range(0,m):
            # state indice is one more in right hand side than in left hand side
            N[0][i] = 1 # for S specie, appears on lhs of every first m reactions
            N[i+1][i] = 1 # diagonal
            M[1][i] = 1 # for I_1 specie, appears on lhs of every first m reactions
            M[i+1][i] += 1 # diagonal + 1 to account for case of M[1][1] where we need 2
        # reactions I_{k-1} -> I_k
        # there are m-1 of them
        for i in range(0,m-1):
            N[1+i][m+i] = 1
            M[1+1+i][m+i] = 1
        
        R = M - N
        params["nbReactions"] = nbReactions
        params["nbStates"] = nbStates
        params['N_matrix'] = N
        params['M_matrix'] = M
        params['R_matrix'] = R
        print('N = \n' , N)
        print('M = \n' , M)
        print('R = \n' , R)
        print("k_plus = ", k_plus)
        print("x_ini = ", x_ini)

    """
    SI model with well mixed hypothesis and intermediate L step
    and erlang waiting time after that
    """
    if (params["model_name"] == "SIAndIntermediateErlangWellMixed"):
        k_SL = params["k_SL"]
        k_LI = params["k_LI"]
        S_ini = params["S_ini"]
        I_ini = params["I_ini"]
        L_ini = params["L_ini"]
        m = params["nb_intermediate_states"]        

        # indices S = 0, I_1 = 1, I_2 = 2... I (=I_m) = m
        # R = indice m +1
        # N_ij = stoechiomatrix coefficients on left hand side of reactions
        # M_ij = stoechiomatrix coefficients on right hand side of reactions
        # i = state, j = reaction
        nbReactions =  m + 1
        # S, L, m * I_{i=1:m}
        nbStates = m + 1 + 1

        # initialize population with known values
        x_ini = np.zeros(nbStates)
        x_ini[0] = S_ini
        x_ini[1] = L_ini
        x_ini[nbStates -1] = I_ini
        params['x_ini'] = x_ini
        
        k_plus = np.zeros(nbReactions)
        # identical results to SIerlangWellMixed
        k_plus[0] = k_SL
        # k_plus[1] = k_LI
        k_plus[1:nbReactions] = k_LI * m
        # k_plus[2:nbReactions] = k_SL * m * S_ini
        params['k_plus'] = k_plus

        N = np.zeros([nbStates, nbReactions])
        M = np.zeros([nbStates, nbReactions])

        # first reaction: S + I_m -> I_m + L (and all intermediate reactions)
        N[0][0] = 1 ; M[1][0] = 1
        N[nbStates-1][0] = 1 ; M[nbStates-1][0] += 1 # add 1 in case m = 1 and we need to have 2 on the right hand side
        
        # reactions L -> I_1
        # and all other I_{k-1} -> I_k        
        for i in range(1,nbReactions):
            # state indice is one more in right hand side than in left hand side
            N[i][i] = 1 
            M[i+1][i] = 1    

        R = M - N
        params["nbReactions"] = nbReactions
        params["nbStates"] = nbStates
        params['N_matrix'] = N
        params['M_matrix'] = M
        params['R_matrix'] = R
        print('N = \n' , N)
        print('M = \n' , M)
        print('R = \n' , R)
        print("k_plus = ", k_plus) 
        print("x_ini = ", x_ini)             

    """
    SI model + interemdiate state with well mixed hypothesis and erlang waiting time
    for the first reaction reaction and we modify it so that 
    that there is no I_m trigger. [USEFUL TO COMPUTE WAITING TIMES]
    S + I   -> I + I_1 (with rate m * k_{SI})
    I_1     -> I_2     (with rate m * k_{SI})
    ...
    I_{m-1} -> I_{m-1}       (with rate m * k_{SI})    
    where m is the number of intermediate states
    """
    if (params["model_name"] == "SIAndIntermediateErlangWellMixedModified"):
        
        k_SL = params["k_SL"]
        k_LI = params["k_LI"]
        S_ini = params["S_ini"]
        I_ini = params["I_ini"]
        L_ini = params["L_ini"]
        m = params["nb_intermediate_states"]        

        # indices S = 0, I_1 = 1, I_2 = 2... I (=I_m) = m
        # R = indice m +1
        # N_ij = stoechiomatrix coefficients on left hand side of reactions
        # M_ij = stoechiomatrix coefficients on right hand side of reactions
        # i = state, j = reaction
        nbReactions =  m + 1
        # S, L, m * I_{i=1:m}
        nbStates = m + 1 + 1

        # initialize population with known values
        x_ini = np.zeros(nbStates)
        x_ini[0] = S_ini
        x_ini[1] = L_ini
        x_ini[nbStates -1] = I_ini
        params['x_ini'] = x_ini
        
        k_plus = np.zeros(nbReactions)
        # identical results to SIerlangWellMixed
        k_plus[0] = k_SL
        k_plus[1] = k_LI
        # k_plus[1:nbReactions] = k_LI * m
        k_plus[2:nbReactions] = k_SL * m * S_ini
        params['k_plus'] = k_plus

        N = np.zeros([nbStates, nbReactions])
        M = np.zeros([nbStates, nbReactions])

        # first reaction: S + I_m -> I_m + L (and all intermediate reactions)
        N[0][0] = 1 ; M[1][0] = 1
        N[nbStates-1][0] = 1 ; M[nbStates-1][0] += 1 # add 1 in case m = 1 and we need to have 2 on the right hand side
        
        # reactions L -> I_1
        # and all other I_{k-1} -> I_k        
        for i in range(1,nbReactions):
            # state indice is one more in right hand side than in left hand side
            N[i][i] = 1 
            M[i+1][i] = 1    

        R = M - N
        params["nbReactions"] = nbReactions
        params["nbStates"] = nbStates
        params['N_matrix'] = N
        params['M_matrix'] = M
        params['R_matrix'] = R
        print('N = \n' , N)
        print('M = \n' , M)
        print('R = \n' , R)
        print("k_plus = ", k_plus)
        print("x_ini = ", x_ini)

    """
    SI model with linear reaction to infect S for the first reaction 
    and erlang wiaitng time 
    S   -> I_1
    I_1 -> I_2     (with rate m * k_{SI})
    ...
    I_{m-1} -> I       (with rate m * k_{SI})    
    where m is the number of intermediate states
    """
    if (params["model_name"] == "SILinearErlangWellMixed"):
        
        k_SI = params["k_SI"]
        S_ini = params["S_ini"]
        I_ini = params["I_ini"]
        m = params["nb_intermediate_states"]

        Lambda_SI = k_SI * m # recovery after 1 event        

        # indices S = 0, I_1 = 1, I_2 = 2... I (=I_m) = m
        # R = indice m +1
        # N_ij = stoechiomatrix coefficients on left hand side of reactions
        # M_ij = stoechiomatrix coefficients on right hand side of reactions
        # i = state, j = reaction
        nbReactions =  m
        # S, m * I_{i=1:m}, 1*R        
        nbStates = m + 1

        # initialize population with known values
        x_ini = np.zeros(nbStates)
        x_ini[0] = S_ini
        x_ini[m] = I_ini
        params['x_ini'] = x_ini

        # first m reactions  = Lambda_SI = k_SI * m        
        k_plus = np.zeros(nbReactions)
        k_plus[0] = k_SI * (m) * S_ini
        k_plus[1:m] = k_SI * (m) * S_ini
        params['k_plus'] = k_plus

        N = np.zeros([nbStates, nbReactions])
        M = np.zeros([nbStates, nbReactions])
        # first reaction: S-> I_1
        N[0][0] = 1 ; M[1][0] = 1        
        # reactions I_1 -> I_2
        #           I_2 -> I_3
        #               ...
        #       I_{m-1} -> I_m = I
        # i ranges from 1 to m-1
        for i in range(1,m):
            # state indice is one more in right hand side than in left hand side
            N[i][i] = 1 #diagonal matrix
            M[i+1][i] = 1
        R = M - N     
             
        params["nbReactions"] = nbReactions
        params["nbStates"] = nbStates
        params['N_matrix'] = N
        params['M_matrix'] = M
        params['R_matrix'] = R
        print('N = \n' , N)
        print('M = \n' , M)
        print('R = \n' , R)
        print("k_plus = ", k_plus)    

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
            # take 1% of total population and round to floor to make sure
            # we work with integer population numbers
            dic_simu_total_pop[n] = dic_simu_total_pop[n]
            # get fraction of confirmed cases
            I = dic_simu_data_covid[n][params["time_start_day"]]['Confirmed']
            # get fraction of deaths
            D = dic_simu_data_covid[n][params["time_start_day"]]['Deaths']
            R = dic_simu_data_covid[n][params["time_start_day"]]['Recovered']
            S = dic_simu_total_pop[n] - I - D - R

            # adjust k_SI to each node population
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
        # get data for entire france
        for n in G.nodes():
            # take 1% of total population and round to floor to make sure
            # we work with integer population numbers
            dic_simu_total_pop[n] = dic_simu_total_pop[n]
            # get fraction of confirmed cases
            I += dic_simu_data_covid[n][params["time_start_day"]]['Confirmed']
            # get fraction of deaths
            D += dic_simu_data_covid[n][params["time_start_day"]]['Deaths']
            R += dic_simu_data_covid[n][params["time_start_day"]]['Recovered']
            S += dic_simu_total_pop[n] - I - D - R

        # get N = total population
        N = S + I + R + D

        # adjust k_SI to each node population
        k_SI = (k_IR + k_ID) * R_0 / N
        k_plus = np.array([k_SI, k_IR, k_ID])
        x_ini = np.array([S, I, R, D])
        params['k_plus'] = k_plus
        params['x_ini'] = x_ini

    """
    model from Marc Lavielle: http://webpopix.org/covidix19.html
    Model to test influence of exposed state to check we computed the right R0
    SIRLD model: S: susceptible I: infectious, R: recovered, L: isolated, D: Deaths
    I -> L means once individuals are infected they isolate
    I -> R means that some individuals are infected and do not isolate 
    but will not die however (because of low viral charge)
    L -> D mean that once people are isolated they can die
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
        S = 0  ; I = 0 ; R_pop = 0 ; D = 0; N_pop = 0;
        # get data for entire france
        for n in G.nodes():
            D += dic_simu_data_covid[n][params["time_start_day"]]['Deaths']
            # R_pop += dic_simu_data_covid[n][params["time_start_day"]]['Recovered']
            N_pop += dic_simu_total_pop[n]
        
        params['total_population'] = N_pop # sueful as some point for rescale and stuff        
        I = params['I_ini']
        # this is not clear how intial number of recovered is related to the real
        # fraction of recovered people. Make it a fit parameter for now
        # R_pop =  np.floor((1/params['alpha']) * R_pop) # real "recovered people"
        R =  params['R_ini']
        L = params['L_ini'] # people who are isolated and will die or recover
        S = N_pop - I - D - R - L
        # adjust k_SI to each node population
        k_plus = np.array([k_SI / N_pop, k_IR, k_IL, k_LD])
        x_ini = np.array([S, I, R, L, D])
        # print("k_plus = ", k_plus)
        params['k_plus'] = k_plus
        params['x_ini'] = x_ini

    """
    model from Marc Lavielle: http://webpopix.org/covidix19.html
    we modify it so that once pepole are isolated S -> 
    Model to test influence of exposed state to check we computed the right R0
    SIRLD model: S: susceptible I: infectious, R: recovered, L: isolated, D: Deaths
    I -> L means once individuals are infected they isolate
    I -> R means that some individuals are infected and do not isolate 
    but will not die however (because of low viral charge)
    L -> D mean that once people are isolated they can die
    """
    if (params["model_name"] == "RealSIRLDmodifiedWellMixed"):
        # load initial state, and reactions values
        k_SI = params["k_SI"]  # A = 0
        k_IR = params["k_IR"]  # A = 2
        k_IL = params["k_IL"]  # A = 3
        k_LR = params["k_LR"]  # A = 4
        k_LD = params["k_LD"]  # A = 5
        
        # indices S = 0, I = 1, R = 2
        # N_ij = stoechiomatrix coefficients on left hand side of reactions
        # M_ij = stoechiomatrix coefficients on right hand side of reactions
        # i = state, j = reaction
        nbReactions = 5
        nbStates = 5
        N = np.zeros([nbStates, nbReactions])
        M = np.zeros([nbStates, nbReactions])

        # 1st reaction: S + I -> 2I
        N[0][0] = 1 ; M[0][0] = 0
        N[1][0] = 1 ; M[1][0] = 1
        N[2][0] = 0 ; M[2][0] = 0
        N[3][0] = 0 ; M[3][0] = 0
        N[4][0] = 0 ; M[4][0] = 0        
        #2nd reaction: I -> R
        N[0][1] = 0 ; M[0][1] = 0
        N[1][1] = 1 ; M[1][1] = 0
        N[2][1] = 0 ; M[2][1] = 1
        N[3][1] = 0 ; M[3][1] = 0
        N[4][1] = 0 ; M[4][1] = 0
        #3rd reaction: I -> L
        N[0][2] = 0 ; M[0][2] = 0
        N[1][2] = 1 ; M[1][2] = 0
        N[2][2] = 0 ; M[2][2] = 0
        N[3][2] = 0 ; M[3][2] = 1
        N[4][2] = 0 ; M[4][2] = 0
        # 4th reaction: L -> R
        N[0][3] = 0 ; M[0][3] = 0
        N[1][3] = 0 ; M[1][3] = 0
        N[2][3] = 0 ; M[2][3] = 1
        N[3][3] = 1 ; M[3][3] = 0
        N[4][3] = 0 ; M[4][3] = 0
        # 5th reaction: L -> D
        N[0][4] = 0 ; M[0][4] = 0
        N[1][4] = 0 ; M[1][4] = 0
        N[2][4] = 0 ; M[2][4] = 0
        N[3][4] = 1 ; M[3][4] = 0
        N[4][4] = 0 ; M[4][4] = 1    
        
        R = M - N
        
        params["nbReactions"] = nbReactions
        params["nbStates"] = nbStates
        params['N_matrix'] = N
        params['M_matrix'] = M
        params['R_matrix'] = R
        
        # print("N_matrix = \n" , N)
        # print("M_matrix = \n" , M)
        # print("R_matrix = \n" , R)
        
        dic_data_covid = nx.get_node_attributes(G, name='data_covid')
        dic_total_pop  = nx.get_node_attributes(G, name='total_pop')
        dic_simu_total_pop = copy.deepcopy(dic_total_pop)
        dic_simu_data_covid = copy.deepcopy(dic_data_covid)
        S = 0 ;  I = 0 ; R_pop = 0 ; D = 0; N_pop = 0;
        # get data for entire france
        for n in G.nodes():
            D += dic_simu_data_covid[n][params["time_start_day"]]['Deaths']
            # R_pop += dic_simu_data_covid[n][params["time_start_day"]]['Recovered']
            N_pop += dic_simu_total_pop[n]
        
        params['total_population'] = N_pop # sueful as some point for rescale and stuff

        I = params['I_ini']
        # this is not clear how intial number of recovered is related to the real
        # fraction of recovered people. Make it a fit parameter for now
        # R_pop =  np.floor((1/params['alpha']) * R_pop) # real "recovered people"
        R =  params['R_ini']
        L = params['L_ini'] # people who are isolated and will die or recover
        S = N_pop - I - D - R - L
        # adjust k_SI to each node population
        k_plus = np.array([k_SI / N_pop, k_IR, k_IL, k_LR, k_LD])
        x_ini = np.array([S, I, R, L, D])
        # print("k_plus = ", k_plus)
        params['k_plus'] = k_plus
        params['x_ini'] = x_ini


    """
    SEIRLD model: S: susceptible, E: exposed (was in contact with infectious 
    but not yet infectious, I: infectious, R: recovered, L: isolated, D: Deaths
    I -> L means once individuals are infected they isolate
    I -> R means that some individuals are infected and do not isolate 
    but will not die however (because of low viral charge)
    L -> R and L -> mean that once people are isolated they can either die or recover
    """
    if (params["model_name"] == "RealSEIRLDWellMixed"):
        # load initial state, and reactions values
        k_SE = params["k_SE"]  # A = 0
        k_EI = params["k_EI"]  # A = 1
        k_IR = params["k_IR"]  # A = 2
        k_IL = params["k_IL"]  # A = 3
        k_LR = params["k_LR"]  # A = 4
        k_LD = params["k_LD"]  # A = 5
        
        # indices S = 0, I = 1, R = 2
        # N_ij = stoechiomatrix coefficients on left hand side of reactions
        # M_ij = stoechiomatrix coefficients on right hand side of reactions
        # i = state, j = reaction
        nbReactions = 6
        nbStates = 6
        N = np.zeros([nbStates, nbReactions])
        M = np.zeros([nbStates, nbReactions])

        # 1st reaction: S + I -> I + E
        N[0][0] = 1 ; M[0][0] = 0
        N[1][0] = 0 ; M[1][0] = 1
        N[2][0] = 1 ; M[2][0] = 1
        N[3][0] = 0 ; M[3][0] = 0
        N[4][0] = 0 ; M[4][0] = 0
        N[5][0] = 0 ; M[5][0] = 0
        # 2nd reaction: E -> I
        N[0][1] = 0 ; M[0][1] = 0
        N[1][1] = 1 ; M[1][1] = 0
        N[2][1] = 0 ; M[2][1] = 1
        N[3][1] = 0 ; M[3][1] = 0
        N[4][1] = 0 ; M[4][1] = 0
        N[5][1] = 0 ; M[5][1] = 0
        # 3rd reaction: I -> R
        N[0][2] = 0 ; M[0][2] = 0
        N[1][2] = 0 ; M[1][2] = 0
        N[2][2] = 1 ; M[2][2] = 0
        N[3][2] = 0 ; M[3][2] = 1
        N[4][2] = 0 ; M[4][2] = 0
        N[5][2] = 0 ; M[5][2] = 0
        # 4th reaction: I -> L
        N[0][3] = 0 ; M[0][3] = 0
        N[1][3] = 0 ; M[1][3] = 0
        N[2][3] = 1 ; M[2][3] = 0
        N[3][3] = 0 ; M[3][3] = 0
        N[4][3] = 0 ; M[4][3] = 1
        N[5][3] = 0 ; M[5][3] = 0
        # 5th reaction: L -> R
        N[0][4] = 0 ; M[0][4] = 0
        N[1][4] = 0 ; M[1][4] = 0
        N[2][4] = 0 ; M[2][4] = 0
        N[3][4] = 0 ; M[3][4] = 1
        N[4][4] = 1 ; M[4][4] = 0
        N[5][4] = 0 ; M[5][4] = 0
        # 6th reaction: L -> D
        N[0][5] = 0 ; M[0][5] = 0
        N[1][5] = 0 ; M[1][5] = 0
        N[2][5] = 0 ; M[2][5] = 0
        N[3][5] = 0 ; M[3][5] = 0
        N[4][5] = 1 ; M[4][5] = 0
        N[5][5] = 0 ; M[5][5] = 1    
        
        R = M - N
        
        params["nbReactions"] = nbReactions
        params["nbStates"] = nbStates
        params['N_matrix'] = N
        params['M_matrix'] = M
        params['R_matrix'] = R
        
        # print("N_matrix = \n" , N)
        # print("M_matrix = \n" , M)
        # print("R_matrix = \n" , R)
        
        dic_data_covid = nx.get_node_attributes(G, name='data_covid')
        dic_total_pop  = nx.get_node_attributes(G, name='total_pop')
        dic_simu_total_pop = copy.deepcopy(dic_total_pop)
        dic_simu_data_covid = copy.deepcopy(dic_data_covid)
        S = 0 ; E = 0 ; I = 0 ; R_pop = 0 ; D = 0; N_pop = 0;
        # get data for entire france
        for n in G.nodes():
            D += dic_simu_data_covid[n][params["time_start_day"]]['Deaths']
            # R_pop += dic_simu_data_covid[n][params["time_start_day"]]['Recovered']
            N_pop += dic_simu_total_pop[n]
        
        params['total_population'] = N_pop # sueful as some point for rescale and stuff
        E = params['E_ini']
        I = params['I_ini']
        # this is not clear how intial number of recovered is related to the real
        # fraction of recovered people. Make it a fit parameter for now
        # R_pop =  np.floor((1/params['alpha']) * R_pop) # real "recovered people"
        R =  params['R_ini']
        L = params['L_ini'] # people who are isolated and will die or recover
        S = N_pop -E - I - D - R - L
        # adjust k_SI to each node population
        k_plus = np.array([k_SE / N_pop, k_EI, k_IR, k_IL, k_LR, k_LD])
        x_ini = np.array([S, E, I, R, L, D])
        # print("k_plus = ", k_plus)
        params['k_plus'] = k_plus
        params['x_ini'] = x_ini

    """
    SEIRLD model on real data (pop, number of cases) and on graph
    we test initialization on grand est and paris for custom plot
    """
    if (params["model_name"] == "RealSEIRLDOnGraphInitTEST"):
        # load initial state, and reactions values
        k_SE = params["k_SE"]  # A = 0
        k_EI = params["k_EI"]  # A = 1
        k_IR = params["k_IR"]  # A = 2
        k_IL = params["k_IL"]  # A = 3
        k_LR = params["k_LR"]  # A = 4
        k_LD = params["k_LD"]  # A = 5
        
        # indices S = 0, I = 1, R = 2
        # N_ij = stoechiomatrix coefficients on left hand side of reactions
        # M_ij = stoechiomatrix coefficients on right hand side of reactions
        # i = state, j = reaction
        nbReactions = 6
        nbStates = 6
        N = np.zeros([nbStates, nbReactions])
        M = np.zeros([nbStates, nbReactions])

        # 1st reaction: S + I -> I + E
        N[0][0] = 1 ; M[0][0] = 0
        N[1][0] = 0 ; M[1][0] = 1
        N[2][0] = 1 ; M[2][0] = 1
        N[3][0] = 0 ; M[3][0] = 0
        N[4][0] = 0 ; M[4][0] = 0
        N[5][0] = 0 ; M[5][0] = 0
        # 2nd reaction: E -> I
        N[0][1] = 0 ; M[0][1] = 0
        N[1][1] = 1 ; M[1][1] = 0
        N[2][1] = 0 ; M[2][1] = 1
        N[3][1] = 0 ; M[3][1] = 0
        N[4][1] = 0 ; M[4][1] = 0
        N[5][1] = 0 ; M[5][1] = 0
        # 3rd reaction: I -> R
        N[0][2] = 0 ; M[0][2] = 0
        N[1][2] = 0 ; M[1][2] = 0
        N[2][2] = 1 ; M[2][2] = 0
        N[3][2] = 0 ; M[3][2] = 1
        N[4][2] = 0 ; M[4][2] = 0
        N[5][2] = 0 ; M[5][2] = 0
        # 4th reaction: I -> L
        N[0][3] = 0 ; M[0][3] = 0
        N[1][3] = 0 ; M[1][3] = 0
        N[2][3] = 1 ; M[2][3] = 0
        N[3][3] = 0 ; M[3][3] = 0
        N[4][3] = 0 ; M[4][3] = 1
        N[5][3] = 0 ; M[5][3] = 0
        # 5th reaction: L -> R
        N[0][4] = 0 ; M[0][4] = 0
        N[1][4] = 0 ; M[1][4] = 0
        N[2][4] = 0 ; M[2][4] = 0
        N[3][4] = 0 ; M[3][4] = 1
        N[4][4] = 1 ; M[4][4] = 0
        N[5][4] = 0 ; M[5][4] = 0
        # 6th reaction: L -> D
        N[0][5] = 0 ; M[0][5] = 0
        N[1][5] = 0 ; M[1][5] = 0
        N[2][5] = 0 ; M[2][5] = 0
        N[3][5] = 0 ; M[3][5] = 0
        N[4][5] = 1 ; M[4][5] = 0
        N[5][5] = 0 ; M[5][5] = 1    
        
        R = M - N
        
        params["nbReactions"] = nbReactions
        params["nbStates"] = nbStates
        params['N_matrix'] = N
        params['M_matrix'] = M
        params['R_matrix'] = R
        
        print("N_matrix = \n" , N)
        print("M_matrix = \n" , M)
        print("R_matrix = \n" , R)
        
        # get dictionnary with data
        dic_data_covid = nx.get_node_attributes(G, name='data_covid')
        dic_total_pop  = nx.get_node_attributes(G, name='total_pop')
        dic_simu_total_pop = copy.deepcopy(dic_total_pop)
        dic_simu_data_covid = copy.deepcopy(dic_data_covid)

        dic_x_ini = {}
        dic_x = {}
        dic_kinetic = {}
        N_pop_france = 0
        for n in G.nodes(): # loop over all compartments                        
            N_pop = dic_simu_total_pop[n]            
            N_pop_france += N_pop
        print("N_pop_france = " , N_pop_france)
        
        density_france = 119.1
        for n in G.nodes(): # loop over all compartments
            D     = dic_simu_data_covid[n][params["time_start_day"]]['Deaths']
            R_pop = dic_simu_data_covid[n][params["time_start_day"]]['Recovered']                        
            N_pop_local = dic_simu_total_pop[n]
            R_pop =  np.floor((1/params['alpha']) * R_pop) # real "recovered people"            
            if (n == 7): # Ile de France
                E = params['E_ini']
                I = params['I_ini']
                L = params['L_ini'] # people who are isolated and will die or recover
                density_local = 1001 #people per km^2
            elif (n == 5) : # Grand-Est
                E = params['E_ini']
                I = params['I_ini']
                L = params['L_ini'] # people who are isolated and will die or recover
                density_local = 97. # people per km^2
            else:
                E = 0
                I = 0                
                L = 0
                # for now: same density everywhere
                density_local = density_france
            
            
            # the constant k_SI given here is \tilde{k_SI_1} = k_SI^1 * N_1
            # where indice 1 refers to that we fitted on France entirely
            # for any other region k_SI^2 = k_SI^1 *(N1/N2) * (n2/n1)
            # the code needs k_SI^2 so we compute 
            # k_SI^2 = \tilde{k_SI^1] / N2 * (n2/n1)
            # TODO: write into notes                
            rescaled_k_SE = k_SE * (density_local / density_france) * (1./N_pop_local)
            
            S = N_pop_local - E - I - D - R_pop - L
            # k_SE is involved in a nonlinear reaction so its value will 
            # be divided by N = total population in the region
            k_plus = np.array([rescaled_k_SE, k_EI, k_IR, k_IL, k_LR, k_LD])
            
            x_ini  = np.array([S, E, I, R_pop, L, D])
            x      = np.array([S, E, I, R_pop, L, D])
            
            dic_x_ini[n] = x_ini
            dic_x[n] = x
            dic_kinetic[n] = k_plus
            
        nx.set_node_attributes(G, name='x_ini', values=dic_x_ini)
        nx.set_node_attributes(G, name='x', values=dic_x)
        nx.set_node_attributes(G, name='k_plus', values=dic_kinetic)
            
        # set edges attributes here
        dic_hopping_rates = {}
        for e in G.edges():
            dic_hopping_rates[e] = 0.0001 # constant for now because TEST
        nx.set_edge_attributes(G, name='hopping_rates', values=dic_hopping_rates)
    
    """
    simple SIR model where we track down every individual
    S + I -> 2I (with rate k_SI)
    I     -> R  (with rate k_IR)
    """
    if (params["model_name"] == "SIRWellMixedAgentBased"):

        k_SI = params["k_SI"]  # A = 1
        k_IR = params["k_IR"]  # A = 2
        S_ini = params["S_ini"]
        I_ini = params["I_ini"]
        R_ini = params["R_ini"]
        
        # declare indices of every state
        x_ini = [[], [], []]
        for i in range(0, S_ini):
            x_ini[0].append(i)
        for i in range(S_ini, S_ini+I_ini):
            x_ini[1].append(i)
        for i in range(S_ini+I_ini, S_ini+I_ini+R_ini):
            x_ini[2].append(i)
        
        params['x_ini'] = copy.deepcopy(x_ini)
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
        
        print("N_matrix = \n" , N)
        print("M_matrix = \n" , M)
        print("R_matrix = \n" , R)
        print('ENDINIT\n\n')

    """
    SI model with well mixed hypothesis and erlang waiting time
    for the first reaction and agent based tracking
    S + I   -> I + I_1 (with rate m * k_{SI})
    I_1     -> I_2     (with rate m * k_{SI})
    ...
    I_{m-1} -> I       (with rate m * k_{SI})    
    where m is the number of intermediate states
    """
    if (params["model_name"] == "SIerlangWellMixedAgentBased"):
        
        k_SI = params["k_SI"]
        S_ini = params["S_ini"]
        I_ini = params["I_ini"]
        m = params["nb_intermediate_states"]

        Lambda_SI = k_SI * m # recovery after 1 event        

        # indices S = 0, I_1 = 1, I_2 = 2... I (=I_m) = m
        # R = indice m +1
        # N_ij = stoechiomatrix coefficients on left hand side of reactions
        # M_ij = stoechiomatrix coefficients on right hand side of reactions
        # i = state, j = reaction
        nbReactions =  m
        # S, m * I_{i=1:m}, 1*R        
        nbStates = m + 1

        # initialize population with known values
        # declare indices of every state
        # intialize list
        x_ini = [[]] # declare list for S
        for k in range(0,m):
            x_ini.append([]) # declare list for all m other states (infected)
        for i in range(0, S_ini):
            x_ini[0].append(i)
        for i in range(S_ini, S_ini+I_ini):
            x_ini[nbStates-1].append(i)        
        
        params['x_ini'] = copy.deepcopy(x_ini)

        # first m reactions  = Lambda_SI
        # second m reactions = Lambda_IR
        k_plus = np.zeros(nbReactions)
        k_plus[0] = k_SI 
        k_plus[1:m] = k_SI *(m-1) * S_ini
        # k_plus[0] = k_SI   * (m)
        # k_plus[1:m] = k_SI * (m)
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

        R = M - N
        params["nbReactions"] = nbReactions
        params["nbStates"] = nbStates
        params['N_matrix'] = N
        params['M_matrix'] = M
        params['R_matrix'] = R
        print('N = \n' , N)
        print('M = \n' , M)
        print('R = \n' , R)
        print("k_plus = ", k_plus) 
        print("x_ini = ", x_ini) 

    """
    SI model with well mixed hypothesis and intermediate L step
    and erlang waiting time after that
    agent based approac
    """
    if (params["model_name"] == "SIAndIntermediateErlangWellMixedAgentBased"):
        k_SL = params["k_SL"]
        k_LI = params["k_LI"]
        S_ini = params["S_ini"]
        I_ini = params["I_ini"]
        L_ini = params["L_ini"]
        m = params["nb_intermediate_states"]        

        # indices S = 0, I_1 = 1, I_2 = 2... I (=I_m) = m
        # R = indice m +1
        # N_ij = stoechiomatrix coefficients on left hand side of reactions
        # M_ij = stoechiomatrix coefficients on right hand side of reactions
        # i = state, j = reaction
        nbReactions =  m + 1
        # S, L, m * I_{i=1:m}
        nbStates = m + 1 + 1

        # initialize population with known values
        # declare indices of every state
        # intialize list
        x_ini = [[], []] # declare list for S and L
        for k in range(0,m):
            x_ini.append([]) # declare list for all m other states (infected)
        for i in range(0, S_ini):
            x_ini[0].append(i)
        # add L initial states
        for i in range(S_ini, S_ini+L_ini):
            x_ini[1].append(i)        
        # add I initial states
        for i in range(S_ini+L_ini, S_ini+L_ini+I_ini):
            x_ini[nbStates-1].append(i)
            
        params['x_ini'] = copy.deepcopy(x_ini)
        
        k_plus = np.zeros(nbReactions)
        # identical results to SIerlangWellMixed
        k_plus[0] = k_SL
        # k_plus[1] = k_LI
        k_plus[1:nbReactions] = k_LI * m
        # k_plus[2:nbReactions] = k_SL * m * S_ini
        params['k_plus'] = k_plus

        N = np.zeros([nbStates, nbReactions])
        M = np.zeros([nbStates, nbReactions])

        # first reaction: S + I_m -> I_m + L (and all intermediate reactions)
        N[0][0] = 1 ; M[1][0] = 1
        N[nbStates-1][0] = 1 ; M[nbStates-1][0] += 1 # add 1 in case m = 1 and we need to have 2 on the right hand side
        
        # reactions L -> I_1
        # and all other I_{k-1} -> I_k        
        for i in range(1,nbReactions):
            # state indice is one more in right hand side than in left hand side
            N[i][i] = 1 
            M[i+1][i] = 1    

        R = M - N
        params["nbReactions"] = nbReactions
        params["nbStates"] = nbStates
        params['N_matrix'] = N
        params['M_matrix'] = M
        params['R_matrix'] = R
        print('N = \n' , N)
        print('M = \n' , M)
        print('R = \n' , R)
        print("k_plus = ", k_plus) 
        print("x_ini = ", x_ini)             

    """
    SI model with well mixed hypothesis and erlang waiting time
    we remove individuals that reach the state I right after so that 
    the number of individual in I stazte is constant    
    agent based approach
    """
    if (params["model_name"] == "SIerlangModifiedWellMixedAgentBased"):
        k_SI = params["k_SI"]
        S_ini = params["S_ini"]
        I_ini = params["I_ini"]
        m = params["nb_intermediate_states"]

        Lambda_SI = k_SI * m # recovery after 1 event        

        # indices S = 0, I_1 = 1, I_2 = 2... I (=I_m) = m
        # R = indice m +1
        # N_ij = stoechiomatrix coefficients on left hand side of reactions
        # M_ij = stoechiomatrix coefficients on right hand side of reactions
        # i = state, j = reaction
        nbReactions =  m
        # S, m * I_{i=1:m}, 1*R        
        nbStates = m + 1

        # initialize population with known values
        # declare indices of every state
        # intialize list
        x_ini = [[]] # declare list for S
        for k in range(0,m):
            x_ini.append([]) # declare list for all m other states (infected)
        for i in range(0, S_ini):
            x_ini[0].append(i)
        for i in range(S_ini, S_ini+I_ini):
            # x_ini[1].append(i)
            x_ini[nbStates-1].append(i)        
        
        params['x_ini'] = copy.deepcopy(x_ini)

        # first m reactions  = Lambda_SI
        # second m reactions = Lambda_IR
        k_plus = np.zeros(nbReactions)
        k_plus[0] = k_SI * m 
        k_plus[1:m] = k_SI *(m)        
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

        R = M - N
        params["nbReactions"] = nbReactions
        params["nbStates"] = nbStates
        params['N_matrix'] = N
        params['M_matrix'] = M
        params['R_matrix'] = R
        print('N = \n' , N)
        print('M = \n' , M)
        print('R = \n' , R)
        print("k_plus = ", k_plus) 
        print("x_ini = ", x_ini) 
    
    """
    SIR model with well mixed hypothesis and another type of
    erlang waiting time (same as in krylova paper)    
    m is the number of intermediate states
    with agent based approach
    """
    if (params["model_name"] == "SIRerlangKrylovaWellMixedAgentBased"):
        
        k_SI = params["k_SI"]
        k_IR = params["k_IR"]
        S_ini = params["S_ini"]
        I_ini = params["I_ini"]
        R_ini = params["R_ini"]
        m = params["nb_intermediate_states"]        

        # indices S = 0, I_1 = 1, I_2 = 2... I (=I_m) = m
        # R = indice m +1
        # N_ij = stoechiomatrix coefficients on left hand side of reactions
        # M_ij = stoechiomatrix coefficients on right hand side of reactions
        # i = state, j = reaction
        nbReactions =  m + m-1 + 1
        # S, m * I_{i=1:m}, 1*R        
        nbStates = m + 1 + 1

        # initialize population with known values
        # declare indices of every state
        # intialize list
        x_ini = [[], []] # declare list for S and R
        for k in range(0,m):
            x_ini.append([]) # declare list for all m other states (infected)
        for i in range(0, S_ini):
            x_ini[0].append(i)        
        for i in range(S_ini, S_ini+I_ini):
            # x_ini[1].append(i) # initialiy all infected are in I_1
            x_ini[nbStates-2].append(i)        
        
        params['x_ini'] = copy.deepcopy(x_ini)

        
        # first m reactions  = Lambda_SI
        # second m reactions = Lambda_IR
        k_plus = np.zeros(nbReactions)
        k_plus[0:m] = k_SI
        k_plus[m:nbReactions] = k_IR * m
        # k_plus[0:m] = k_SI * m
        # k_plus[m:nbReactions] = k_SI * m
        params['k_plus'] = k_plus

        N = np.zeros([nbStates, nbReactions])
        M = np.zeros([nbStates, nbReactions])
        
        # reactions S + I_[k} -> I_1 + I_k
        # there are m of them
        for i in range(0,m):
            # state indice is one more in right hand side than in left hand side
            N[0][i] = 1 # for S specie, appears on lhs of every first m reactions
            N[i+1][i] = 1 # diagonal
            M[1][i] = 1 # for I_1 specie, appears on lhs of every first m reactions
            M[i+1][i] += 1 # diagonal + 1 to account for case of M[1][1] where we need 2
        # reactions I_{k-1} -> I_k
        # there are m-1 of them
        for i in range(0,m-1):
            N[1+i][m+i] = 1
            M[1+1+i][m+i] = 1
        
        # last reaction: I_m --> R
        N[nbStates-2][nbReactions-1] = 1
        M[nbStates-1][nbReactions-1] = 1
        
        R = M - N
        params["nbReactions"] = nbReactions
        params["nbStates"] = nbStates
        params['N_matrix'] = N
        params['M_matrix'] = M
        params['R_matrix'] = R
        print('N = \n' , N)
        print('M = \n' , M)
        print('R = \n' , R)
        print("k_plus = ", k_plus)
        print("x_ini = ", x_ini)

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

    if (params["model_name"] == "SIerlangWellMixedModified") or (params["model_name"] == "SIAndIntermediateErlangWellMixedModified"):
        waiting_times = []

    x = x_ini
    t = params["time_start"]
    outputs = []  # we will fill this list with data to be plotted
    count_loops = 0
    try: # so we can exit the loop at any time with ctrl + C
        while t < params["time_end"]:
            if ((count_loops % params['save_frequency']) == 0):
                x_save = copy.deepcopy(x)
                outputs.append((t, x_save))
                if ((count_loops % 100000) == 0):            
                    print(t,count_loops, x_save)
                # print(t,count_loops, x_save)
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

            # print("t_forward_plus = " , t_forward_plus)
            # print("x = " , x)
            # total rate which the state can be left with
            W = np.sum(t_forward_plus)
            xi = randZeroOne()
            # print(t,count_loops, x, W, t_forward_plus)
            # xi = np.random.uniform(0,1)
            # if the total rate is 0, nothing more can happen anymore
            if (abs(W) <= 0.0000000001):
                # print("W = " , W)
                # print("x = " , x)
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
            # xi = np.random.uniform(0,1)
            if xi <= next_states[0]:
                id_reaction = 0
            else:
                for j in range(1, nbReactions):
                    if (next_states[j - 1] < xi) and (xi < next_states[j]):
                        id_reaction = j
            
            
            # set new state according to reaction
            
            x = x + R[:, id_reaction]
            # special case of a custom model to test waiting times...
            if (params["model_name"] == "SIerlangWellMixedModified") or (params["model_name"] == "SIAndIntermediateErlangWellMixedModified"):
                # if we triggered the last reaction, meaning I_m was increased by 1
                # we register the time at which it happened and remove it right away
                if (id_reaction == nbReactions-1) :
                    waiting_times.append(t)
                    x[nbStates - 1] -= 1
            # print(count_loops, t, x, t_forward_plus, next_states, id_reaction)
            # print(count_loops, t, x, t_forward_plus, W)

            count_loops += 1        
        
    
    except KeyboardInterrupt:
        pass
    
    print(t,count_loops, x)
    
    if (params["model_name"] == "SIerlangWellMixedModified") or (params["model_name"] == "SIAndIntermediateErlangWellMixedModified"):
        return outputs , waiting_times
    else:
        return outputs
            
def SolveAgentBasedOneCompartment(params):
    """
    solve master equation with one compartment
    and track every individual in a agent based fashion
    if the rate of the reaction for state S e.g. is k then every
    individual S^1, S^2, ...S^n such that 
    S^1 + S^2 + ... + S^n = S (total number of individuals)
    each individual can jump from state S^1 to state I^1    
    """
    nbReactions   = params["nbReactions"]
    nbStates      = params["nbStates"]
    k_plus        = params['k_plus']
    N             = params['N_matrix']
    R             = params['R_matrix']
    # x_ini is a list of lists of size nbStates, each sub-list contains
    # the individuals indices in each state
    x_ini          = params['x_ini']
    indices_states_per_reaction = np.zeros([nbStates, nbReactions])
    for i in range(0,nbStates):
        for j in range(0,nbReactions):
            # we don't know how to deal with such reactions yet
            # example : 2S + I --> 3I
            # how does ont translate forward rate t_0 = k_SI * S * (S-1) * I
            # when dealing with individuals??
            if (abs(R[i,j]) > 1) :
                print('ERROR: the R_matrix given that accounts for the reaction changes is of too high order')
                print('R[i,j] can only be = 1 or -1')
                print('We foud R[' , i , ',' , j ,'] = ' , R[i,j])
            
            # fill matrix that accuonts if a state is in a given reaction
            if (abs(N[i,j]) > 0.00000001):  
                indices_states_per_reaction[i,j] = 1
                
    if (params["model_name"] == "SIerlangModifiedWellMixedAgentBased"):
        waiting_times = []
                
    x = copy.deepcopy(x_ini)
    t = params["time_start"]
    outputs = []  # we will fill this list with data to be plotted
    count_loops = 0
    try: # so we can exit the loop at any time with ctrl + C
        while t < params["time_end"]:
            if ((count_loops % params['save_frequency']) == 0):
                x_save = copy.deepcopy(x)
                outputs.append((t, x_save))
                if ((count_loops % 100000) == 0):            
                    print(t,count_loops, x_save)
                # print(t,count_loops, x_save)                        
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
                            y = y * (len(x[i]) -k) # x[i] is a list with all indices that are in 

                    t_forward_plus[j] = t_forward_plus[j] * y

                # finally multiply the forward reaction rate by the reaction rate k+
                t_forward_plus[j] = t_forward_plus[j] * k_plus[j]

            # print("t_forward_plus = " , t_forward_plus)
            # print("x = " , x)
            # total rate which the state can be left with
            W = np.sum(t_forward_plus)
            xi = randZeroOne()
            # print(t,count_loops, x, W, t_forward_plus)
            # xi = np.random.uniform(0,1)
            # if the total rate is 0, nothing more can happen anymore
            if (abs(W) <= 0.0000000001):
                # print("W = " , W)
                # print("x = " , x)
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
            # xi = np.random.uniform(0,1)
            if xi <= next_states[0]:
                id_reaction = 0
            else:
                for j in range(1, nbReactions):
                    if (next_states[j - 1] < xi) and (xi < next_states[j]):
                        id_reaction = j
            
            # we chose reaction id_reaction but we now need to know to 
            # which sub-reaction (i.e. which individual in particular 
            # moves)
            # compute total number of sub-reactions. e.g. t_0 = k * S * I
            # that would be k * S1 * I1, k * S2 * I1, ..., k * Sn * I1,
            # k * S1 * I2 , ... * k Sn * I2
            nbStatesInReaction = 0
            nbSubReactions = 1
            for kk in range(0, nbStates):
                # get number of combinations that can happen between each pair 
                # of states involve in each reaction
                if (abs(indices_states_per_reaction[kk][id_reaction]) > 0.0000001):
                    nbSubReactions = nbSubReactions * indices_states_per_reaction[kk][id_reaction] * len(x[kk])
                    nbStatesInReaction += indices_states_per_reaction[kk][id_reaction]
            
            # each sub event has the same proba and their sum equals
            # t_forward_plus[id_reaction] 
            # we keep it normalized to match the value of xi which is 
            # between 0 and 1                        
            proba_each_sub_event = t_forward_plus[id_reaction] / nbSubReactions / W
            # print("id_reaction, nbStatesInReaction, nbSubReactions = " , id_reaction, nbStatesInReaction, nbSubReactions)
            if (id_reaction == 0): # first reaction we loop from 0
                min_lbound = 0.
            else: # start from the previous value of the chosen reaction
                min_lbound = next_states[id_reaction -1]
            indice = 0
            max_ubound = next_states[id_reaction]
            lbound = min_lbound
            ubound = min_lbound + proba_each_sub_event    
            done = False
            if (xi < lbound): # lucky, we found the first indice as quickly as possible
                indice = 0 
            else: # we need to iterate between next_states[id_reaction-1] and next_states[id_reaction] 
                while (not(done)):
                    # we found the sub event that 
                    if (lbound < xi) and (xi < ubound):
                        done = True
                    elif (xi > max_ubound):
                        print('Something went wrong we went over the interval where we were supposed to look for the indice of individual to move. Exit here')
                        exit(0)
                    else :
                        # increase bounds by one sub event
                        lbound += proba_each_sub_event
                        ubound += proba_each_sub_event
                        indice += 1                    
            
            # in case where 1 state is involved in reaction id_reaction 
            # the indice we found is the same as the position in the corresponding
            # vector of states
            
            if (nbStatesInReaction == 1):
                # assume sub-reactions are ordered as follow:
                # k * x[state][0], k * x[state][1], ... , k * x[state][n]                
                indice_sub_event_state_to_remove = indice
                for ii in range(0,nbStates):
                    if (R[ii, id_reaction] == -1):
                        state_to_remove = ii
                    elif (R[ii, id_reaction] == 1):
                        state_to_add = ii                                                        
                # get id of individual we will remove
                id_individual_state_to_remove = x[state_to_remove][indice_sub_event_state_to_remove]
                # remove individual from list
                x[state_to_remove].pop(indice_sub_event_state_to_remove)
                # add individual to state to add
                x[state_to_add].append(id_individual_state_to_remove)
                        
            if (nbStatesInReaction == 2):
                # assume sub-reactions are ordered as follow: with 
                # state1 = state_to_remove and state2 = state_to_add
                # example for 3 individuals in state 1 and 2 
                # individuals in state 2
                # k * x[state1][0] * x[state2][0], k * 
                # x[state1][1] * x[state2][0], k * x[state1][2] * 
                # x[state2][0], k * x[state1][0] * x[state2][1], k * 
                # x[state1][1] * x[state2][1], , k * x[state1][2] * 
                # x[state2][1]
                # 0 0, 1 0, 2 0, 0 1, 1 1, 2 1
                # found_first_state = False
                for ii in range(0,nbStates):                  
                    if ((R[ii, id_reaction]) == -1):
                        state_to_remove = ii
                    elif ((R[ii, id_reaction]) == 1):
                        state_to_add = ii
                
                # get indice of the sub event
                indice_sub_event_state_to_remove = indice % len(x[state_to_remove])
                                
                # get id of individual we will remove
                id_individual_state_to_remove = x[state_to_remove][indice_sub_event_state_to_remove]
                # remove individual from list
                x[state_to_remove].pop(indice_sub_event_state_to_remove)
                # add individual to state to add
                x[state_to_add].append(id_individual_state_to_remove)
            
            # in the case of erlang modified we remove individuals 
            # that reached infected state right after they did and 
            # register the corresponding waiting time
            if (params["model_name"] == "SIerlangModifiedWellMixedAgentBased"):
                # if we triggered the last reaction, meaning I_m was increased by 1
                # we register the time at which it happened and remove it right away
                if (id_reaction == nbReactions-1) :
                    waiting_times.append(t)
                    # remove individual that we just added to the new state
                    # print("id_reaction, id_individual_state_to_remove = " , id_reaction, id_individual_state_to_remove)
                    x[state_to_add].remove(id_individual_state_to_remove)
                    
            # print('moved ' , id_individual_state_to_remove , 'from' , state_to_remove, 'to ', state_to_add)
            # print(count_loops, t, x, t_forward_plus, next_states, id_reaction)
            # print(count_loops, t, x, id_reaction)

            count_loops += 1

    except KeyboardInterrupt:
        pass
    
    print(t,count_loops, x)
    
    if (params["model_name"] == "SIerlangModifiedWellMixedAgentBased") :
        return outputs , waiting_times
    else:
        return outputs
         
            
    return outputs

def SolveOnGraph(params, G):
    """
    solve master equation with several compartments (nodes of graph)
    with hopping processes
    """
    nbReactions = params["nbReactions"]
    nbStates    = params["nbStates"]
    N           = params['N_matrix']
    R           = params['R_matrix']

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
                print("t, count_loops, x = " , t, count_loops, x_print)
                # if ((count_loops % 10000) == 0):
                    # print("t, count_loops, x = " , t, count_loops, x_print)
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
            # print("W_inside, W_between, W_total = " , W_inside, W_between, W_total)
                
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
                                    if (i > 0):
                                        print('a hopping happens from someone who IS NOT susceptible i = ' , i)                                    
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
