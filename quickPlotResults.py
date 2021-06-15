#!/usr/bin/env python3

from numpy import *
import matplotlib.pyplot as plt
from plotModels import plot_SILinearErlangWellMixed
import collections
import numpy as np

### -----------------------------------------------------------------###
### -----------------------------------------------------------------###
### --- case SI model with S fixed, erlang with m =2
### -----------------------------------------------------------------###
### -----------------------------------------------------------------###
## _dir='/home/jdichamp/INRIA/covid19_epidemic_modeling/cervoise/outputs/VALIDATIONS/SIerlang_one_compartment/rates_k_SI_and_k_SI_times_m_times_S_ini/'
#_dir='/home/jdichamp/INRIA/covid19_epidemic_modeling/cervoise/outputs/VALIDATIONS/SIerlangFixS_one_compartment/save/'
## fname1=_dir+'outputs_SIerlangWellMixed_m=1.txt'
#fname1=_dir+'outputs_SIerlangFixSWellMixed_m=1.txt'
#s = open(fname1, 'r').read()
#outputs_m1 = list(eval(s))

## params1_name = _dir+'params_SIerlangWellMixed_m=1.txt'
#params1_name = _dir+'params_SIerlangFixSWellMixed_m=1.txt'
#s = open(params1_name, 'r').read()
#params_m1 = dict(eval(s))

#_dir='/home/jdichamp/INRIA/covid19_epidemic_modeling/cervoise/outputs/VALIDATIONS/SIerlangFixS_one_compartment/'
## params10_name = _dir+'params_SIerlangWellMixed_m=10.txt'
#params10_name = _dir+'params_SIerlangFixSWellMixed_m=2.txt'
#s = open(params10_name, 'r').read()
#params_m2 = dict(eval(s))

## fname10=_dir+'outputs_SIerlangWellMixed_m=10.txt'
#fname10=_dir+'outputs_SIerlangFixSWellMixed_m=2.txt'
#s = open(fname10, 'r').read()
#outputs_m2 = list(eval(s))


#nbTimeSteps_m1 = len(outputs_m1)
#time_m1=np.zeros(nbTimeSteps_m1)
#nbStates_m1 = outputs_m1[0][1].shape[0]
#m1 = params_m1["nb_intermediate_states"]
#x_m1 = np.zeros([nbTimeSteps_m1, nbStates_m1])
#for i in range(0,nbTimeSteps_m1):
    #time_m1[i] = outputs_m1[i][0]
    #x_m1[i,:] = outputs_m1[i][1]        

#nbTimeSteps_m2 = len(outputs_m2)
#time_m2=np.zeros(nbTimeSteps_m2)
#nbStates_m2 = outputs_m2[0][1].shape[0]
#m2 = params_m2["nb_intermediate_states"]
#x_m2 = np.zeros([nbTimeSteps_m2, nbStates_m2])
#for i in range(0,nbTimeSteps_m2):
    #time_m2[i] = outputs_m2[i][0]
    #x_m2[i,:] = outputs_m2[i][1]        

#tmp_sum = np.zeros(nbTimeSteps_m2)
#for ii in range(0,nbStates_m2):
    #if (ii ==0):
        #tmp_sum += x_m2[:,ii]
    #if((ii>0) and (ii<m2)):
        #tmp_sum += x_m2[:,ii]
#print(nbTimeSteps_m1, nbTimeSteps_m2)
#print("m = 1, slope I = " , (np.log(x_m1[nbTimeSteps_m1-1,1]) - np.log(x_m1[int(nbTimeSteps_m1/2.),1]))/(time_m1[nbTimeSteps_m1-1] - time_m1[int(nbTimeSteps_m1/2.)]))
#print("m = 2, slope I1 = " , (np.log(x_m2[nbTimeSteps_m2-1,1]) - np.log(x_m2[int(nbTimeSteps_m2/2.),1]))/(time_m2[nbTimeSteps_m2-1] - time_m2[int(nbTimeSteps_m2/2.)]))
#print("m = 2, slope I2 = " , (np.log(x_m2[nbTimeSteps_m2-1,2]) - np.log(x_m2[int(nbTimeSteps_m2/2.),2]))/(time_m2[nbTimeSteps_m2-1] - time_m2[int(nbTimeSteps_m2/2.)]))

#fig1 = plt.figure(figsize=(16, 10)) # in inches!
#plt.plot(time_m1, x_m1[:,0], label='S (m=1)')
#plt.plot(time_m1, x_m1[:,nbStates_m1-1], label='I (m=1)')
#plt.plot(time_m2, x_m2[:,1], label='I1 (m=2)')
#plt.plot(time_m2, x_m2[:,1] + x_m2[:,2], label='I1 + I2 (m=2)')
#plt.plot(time_m2, x_m2[:,2], label='I2 (m=2)')
#plt.yscale('log')
#plt.legend()
#plt.xlabel('Arbitrary time unit')
#plt.ylabel('Population numbers')
#plt.show()

### -----------------------------------------------------------------###
### -----------------------------------------------------------------###
### --- case SI model, erlang with m =2
### -----------------------------------------------------------------###
### -----------------------------------------------------------------###

# _dir='/home/jdichamp/INRIA/covid19_epidemic_modeling/cervoise/outputs/VALIDATIONS/SIerlang_one_compartment/'
# fname1=_dir+'outputs_SIerlangWellMixed_m=1.txt'
# s = open(fname1, 'r').read()
# outputs_m1 = list(eval(s))

# params1_name = _dir+'params_SIerlangWellMixed_m=1.txt'
# s = open(params1_name, 'r').read()
# params_m1 = dict(eval(s))

# _dir='/home/jdichamp/INRIA/covid19_epidemic_modeling/cervoise/outputs/VALIDATIONS/SIerlang_one_compartment/'
# params10_name = _dir+'params_SIerlangWellMixed_m=2.txt'
# s = open(params10_name, 'r').read()
# params_m2 = dict(eval(s))

# fname10=_dir+'outputs_SIerlangWellMixed_m=2.txt'
# s = open(fname10, 'r').read()
# outputs_m2 = list(eval(s))


# nbTimeSteps_m1 = len(outputs_m1)
# time_m1=np.zeros(nbTimeSteps_m1)
# nbStates_m1 = outputs_m1[0][1].shape[0]
# m1 = params_m1["nb_intermediate_states"]
# x_m1 = np.zeros([nbTimeSteps_m1, nbStates_m1])
# for i in range(0,nbTimeSteps_m1):
    # time_m1[i] = outputs_m1[i][0]
    # x_m1[i,:] = outputs_m1[i][1]        

# nbTimeSteps_m2 = len(outputs_m2)
# time_m2=np.zeros(nbTimeSteps_m2)
# nbStates_m2 = outputs_m2[0][1].shape[0]
# m2 = params_m2["nb_intermediate_states"]
# x_m2 = np.zeros([nbTimeSteps_m2, nbStates_m2])
# for i in range(0,nbTimeSteps_m2):
    # time_m2[i] = outputs_m2[i][0]
    # x_m2[i,:] = outputs_m2[i][1]        

# tmp_sum = np.zeros(nbTimeSteps_m2)
# for ii in range(0,nbStates_m2):
    # if (ii ==0):
        # tmp_sum += x_m2[:,ii]
    # if((ii>0) and (ii<m2)):
        # tmp_sum += x_m2[:,ii]

# fig1 = plt.figure(figsize=(16, 10)) # in inches!
# plt.plot(time_m1, x_m1[:,0], label='S (m=1)')
# plt.plot(time_m1, x_m1[:,1], label='I1 (m=1)')
# plt.plot(time_m2, x_m2[:,0], label='S (m=2)')
# plt.plot(time_m2, x_m2[:,1], label='I1 (m=2)')
# plt.plot(time_m2, x_m2[:,2], label='I2 (m=2)')
# plt.plot(time_m2, x_m2[:,0] + x_m2[:,1], label='S + I1 (m=2)')
# plt.yscale('log')
# plt.legend()
# plt.xlabel('Arbitrary time unit')
# plt.ylabel('Population numbers')
# plt.show()



### -----------------------------------------------------------------###
### -----------------------------------------------------------------###
### --- case SI model with L interemdiate compartment, erlang with m =10
### -----------------------------------------------------------------###
### -----------------------------------------------------------------###

_dir='/home/jdichamp/INRIA/covid19_epidemic_modeling/cervoise/outputs/VALIDATIONS/SIerlang_one_compartment/'
fname1=_dir+'outputs_SIAndIntermediateErlangWellMixed_m=1.txt'
s = open(fname1, 'r').read()
outputs_m1 = list(eval(s))

params1_name = _dir+'params_SIAndIntermediateErlangWellMixed_m=1.txt'
s = open(params1_name, 'r').read()
params_m1 = dict(eval(s))

_dir='/home/jdichamp/INRIA/covid19_epidemic_modeling/cervoise/outputs/VALIDATIONS/SIerlang_one_compartment/'
params10_name = _dir+'params_SIAndIntermediateErlangWellMixed_m=10.txt'
s = open(params10_name, 'r').read()
params_m10 = dict(eval(s))

fname10=_dir+'outputs_SIAndIntermediateErlangWellMixed_m=10.txt'
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
m2 = params_m10["nb_intermediate_states"]
x_m10 = np.zeros([nbTimeSteps_m10, nbStates_m10])
for i in range(0,nbTimeSteps_m10):
    time_m10[i] = outputs_m10[i][0]
    x_m10[i,:] = outputs_m10[i][1]        

sum_S_I_i_case_m10_and_intermediate_L = np.zeros(nbTimeSteps_m10)
for ii in range(0,nbStates_m10-1):    
    sum_S_I_i_case_m10_and_intermediate_L += x_m10[:,ii]        

fig1 = plt.figure(figsize=(16, 10)) # in inches!
# plt.plot(time_m1, x_m1[:,0], label='S (m=1)')
plt.plot(time_m1, x_m1[:,0] + x_m1[:,1], label='S + L (m=1)')
plt.plot(time_m1, x_m1[:,2], label='I1 (m=1)')
# plt.plot(time_m10, x_m10[:,11], label='S (m=10, with L)')
plt.plot(time_m10, sum_S_I_i_case_m10_and_intermediate_L, label='S + L + I_i (i=1:m-1) (m=10, with L)')
plt.plot(time_m10, x_m10[:,11], label='I10 (m=10, with L)')
# plt.yscale('log')
plt.legend()
plt.xlabel('Arbitrary time unit')
plt.ylabel('Population numbers')
plt.show()

