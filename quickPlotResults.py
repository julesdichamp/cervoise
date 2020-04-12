#!/usr/bin/env python3

from numpy import *
import matplotlib.pyplot as plt
from plotModels import plot_SILinearErlangWellMixed
import collections
import numpy as np


_dir='/home/jdichamp/INRIA/covid19_epidemic_modeling/cervoise/outputs/VALIDATIONS/SIerlang_one_compartment/rates_k_SI_and_k_SI_times_m_times_S_ini/'
fname1=_dir+'outputs_SIerlangWellMixed_m=1.txt'
s = open(fname1, 'r').read()
outputs_m1 = list(eval(s))

params1_name = _dir+'params_SIerlangWellMixed_m=1.txt'
s = open(params1_name, 'r').read()
params_m1 = dict(eval(s))

params10_name = _dir+'params_SIerlangWellMixed_m=10.txt'
s = open(params10_name, 'r').read()
params_m10 = dict(eval(s))

fname10=_dir+'outputs_SIerlangWellMixed_m=10.txt'
s = open(fname10, 'r').read()
outputs_m10 = list(eval(s))


nbTimeSteps_m1 = len(outputs_m1)
time_m1=np.zeros(nbTimeSteps_m1)
nbStates_m1 = outputs_m1[0][1].shape[0]
m1 = params_m1["nb_intermediate_states"]
x_m1 = np.zeros([nbTimeSteps_m1, nbStates_m1])
for i in range(0,nbTimeSteps_m1):
    time_m1[i] = outputs_m1[i][0]
    x_m1[i,:] = outputs_m1[i][1]        

nbTimeSteps_m10 = len(outputs_m10)
time_m10=np.zeros(nbTimeSteps_m10)
nbStates_m10 = outputs_m10[0][1].shape[0]
m10 = params_m10["nb_intermediate_states"]
x_m10 = np.zeros([nbTimeSteps_m10, nbStates_m10])
for i in range(0,nbTimeSteps_m10):
    time_m10[i] = outputs_m10[i][0]
    x_m10[i,:] = outputs_m10[i][1]        

tmp_sum = np.zeros(nbTimeSteps_m10)
for ii in range(0,nbStates_m10):
    if (ii ==0):
        tmp_sum += x_m10[:,ii]
    if((ii>0) and (ii<m10)):
        tmp_sum += x_m10[:,ii]
    
N = 18180
fig1 = plt.figure(figsize=(16, 10)) # in inches!
plt.plot(time_m1/(np.log(2)**2), x_m1[:,0], label='S (m=1)')
plt.plot(time_m1/(np.log(2)**2), x_m1[:,nbStates_m1-1], label='I (m=1)')
plt.plot(time_m10,tmp_sum,color='gray',label='sum of S + I_i (i=1:m-1) (m=10)')
plt.plot(time_m10,x_m10[:,0],label='S (m=10)')
plt.plot(time_m10, x_m10[:,nbStates_m10-1], label='I (m=10)')
# plt.yscale('log')
plt.legend()
plt.xlabel('Arbitrary time unit')
plt.ylabel('Population numbers')
plt.show()

# plot_SILinearErlangWellMixed(params_m1, outputs_m1)

# plt.plot(num_iter,fitness)
# plt.yscale('log')
# plt.show()
