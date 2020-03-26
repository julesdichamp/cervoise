#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import geopandas as gpd
import pandas as pd
import numpy as np


def SIRWellMixed(params):
	"""
	simple SIR model: 
	Susceptible -> Infectious
	Infectious -> Recovered
	k_SI is rate from Susceptible to Infectious
	k_IR is rate from Infectious to Recovered
	"""
	print('SIRWellMixed')
	
	# load initial state, and reactions values
	k_SI = params['k_SI'] # A = 1
	k_IR = params['k_IR'] # A = 2
	S_ini = params['S_ini']
	I_ini = params['I_ini']
	R_ini = params['R_ini']
		
	x_ini = np.array([S_ini, I_ini, R_ini])	
	k_plus = np.array([k_SI, k_IR])
	# indices S = 0, I = 1, R = 2
	# N_ij = stoechiomatrix coefficients on left hand side of reactions
	# M_ij = stoechiomatrix coefficients on right hand side of reactions
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
			
	# intialize states and time
	x = x_ini	
	t = params['time_start']
	outputs = [] # we will fill this list with data to be plotted
	while (t < params['time_end']):
		if (x[0] == 0): # leave loop if there are no more susceptibles
			break				
		# print(t, x)		
		outputs.append((t,x))				
		# forward reaction rates	
		t_forward_plus = np.zeros(nbReactions)
		for j in range(0,nbReactions):
			t_forward_plus[j] = 1.0
			for i in range(0,nbStates):
				t_forward_plus[j] = t_forward_plus[j] * np.math.factorial(x[i]) / np.math.factorial(x[i] - N[i,j])
			t_forward_plus[j] = t_forward_plus[j] * k_plus[j]	
				
		# total rate which the state can be left with	
		W = np.sum(t_forward_plus) 	
		xi = np.random.uniform(0,1)
		
		# Gillepsie time of next event
		dt = - 1 / W * np.log(1-xi)
		t = t + dt
		
		# sum reactions forward rates to choose reaction
		next_states=np.zeros(nbReactions)
		next_states[0] =  t_forward_plus[0]
		for j in range(1,nbReactions):
			next_states[j] = next_states[j-1] + t_forward_plus[j]
		next_states = next_states/W
		
		# choose reaction 
		xi = np.random.uniform(0,1)				
		if (xi <= next_states[0]):
			id_reaction = 0
		else:
			for j in range(1,nbReactions):				
				if ((next_states[j-1] < xi) and (xi < next_states[j]) ):
					id_reaction = j
		
		# set new state according to reaction
		x = x - R[:,id_reaction]

	return outputs
	
def SIRonRegions(params):
	"""
	TODO
	"""	

def SIRonDepartments(params):
	"""
	TODO
	"""



