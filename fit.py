import lmfit
import numpy as np
from numpy.linalg import norm
import pandas as pd
from models import InitModels, SolveOnGraph, SolveOneCompartment
# from datetime import date
import datetime
import networkx as nx

def FitFuncOnGraph(paramsToFit, params, G): 
    print('FitFuncOnGraph')  
    # add parameters to fit into (that will be varied every call to the fitfunction)
    # to the dictionnary that contains all parameters
    params["R_0"]  = paramsToFit["R_0"]
    params["k_IR"] = paramsToFit["k_IR"]
    params["k_ID"] = paramsToFit["k_ID"]    

    # transform to datetime format
    dates = [datetime.datetime.strptime(ts, "%Y-%m-%d") for ts in params["dates_for_fit"]] 
    time_start_date = datetime.datetime.strptime(params['time_start_day'], "%Y-%m-%d")    
    
    # init model and solve
    [params , G] = InitModels(params, G)
    [outputs , G] = SolveOnGraph(params, G)    
    
    nbDates=len(params["dates_for_fit"])
    nbTimes = len(outputs)
    fitness = np.zeros(3 * nbDates * G.number_of_nodes())
    
    # fill vector with relative increase from one date to the other    
    # dates[0] should NOT be the same than time_start_date otherwise 
    # delta_days[0] = 0 and the first fit data point is the intial values...
    # dates bmust be sorted in params['dates_for_fit']
    delta_days = []
    delta_days.append((dates[0] - time_start_date).days) # transform to number of days
    for i in range(1,nbDates):
        delta_days.append((dates[i] - dates[i-1]).days) # transform to number of days
    
    # time_end must verify time_end = time_start + sum(delta_days)
    if params['time_end'] != params['time_start'] + sum(delta_days):        
        params['time_end'] = params['time_start'] + sum(delta_days)
        print("new time_end = " , params['time_end'])

    
    time_date_target = params['time_start'] + delta_days[0]
    _count = 0
    
    # get all data we need    
    dic_data_covid = nx.get_node_attributes(G, name='data_covid')
    # loop over time steps to get simu and compute fit function
    for i in range(0,nbTimes):
        mytime = outputs[i][0]
        # not sure of what precision to take
        if (abs(mytime - time_date_target) <= 0.01) :
            print('[TARGET FOUND]: mytime, date_target = ' , mytime, time_date_target)
            dic_x = outputs[i][1]
            for n in G.nodes():
                # Infected
                fitness[n + _count + 0] = (dic_data_covid[n][params["dates_for_fit"][_count]]['Confirmed']\
                                        - dic_x[n][1])
                fitness[n + _count + 1] = (dic_data_covid[n][params["dates_for_fit"][_count]]['Recovered']\
                                        - dic_x[n][2])
                fitness[n + _count + 2] = (dic_data_covid[n][params["dates_for_fit"][_count]]['Deaths']\
                                        - dic_x[n][3])
                print(n, dic_data_covid[n][params["dates_for_fit"][_count]]['Confirmed'], dic_data_covid[n][params["dates_for_fit"][_count]]['Recovered'], dic_data_covid[n][params["dates_for_fit"][_count]]['Deaths'], dic_x[n][1], dic_x[n][2], dic_x[n][3])
            # get new date to target
            time_date_target = time_date_target + delta_days[_count]
            _count += 1
            # exist loop if we got all data we need
            if (_count == nbDates):
                break;

    print("fitness value = ", norm(fitness, ord=2) ** 2)
    return norm(fitness, ord=2) ** 2

def FitFuncOneComp(paramsToFit, params, G): 
    print('FitFuncOneComp')  
    # add parameters to fit into (that will be varied every call to the fitfunction)
    # to the dictionnary that contains all parameters
    params["R_0"]  = paramsToFit["R_0"]
    params["k_IR"] = paramsToFit["k_IR"]
    params["k_ID"] = paramsToFit["k_ID"]    

    # transform to datetime format
    dates = [datetime.datetime.strptime(ts, "%Y-%m-%d") for ts in params["dates_for_fit"]] 
    time_start_date = datetime.datetime.strptime(params['time_start_day'], "%Y-%m-%d")        
    nbDates=len(params["dates_for_fit"])    
        
    # fill vector with relative increase from one date to the other    
    # dates[0] should NOT be the same than time_start_date otherwise 
    # delta_days[0] = 0 and the first fit data point is the intial values...
    # dates bmust be sorted in params['dates_for_fit']
    delta_days = []
    delta_days.append((dates[0] - time_start_date).days) # transform to number of days
    for i in range(1,nbDates):
        delta_days.append((dates[i] - dates[i-1]).days) # transform to number of days
    
    # time_end must verify time_end = time_start + sum(delta_days)
    if params['time_end'] != params['time_start'] + sum(delta_days):        
        params['time_end'] = params['time_start'] + sum(delta_days)
        print("new time_end = " , params['time_end'])
    
    # init model and solve
    [params , G] = InitModels(params, G)
    outputs = SolveOneCompartment(params)
    nbTimes = len(outputs)
    
    fitness = np.zeros(3 * nbDates)
    time_date_target = params['time_start'] + delta_days[0]
    _count = 0
    # get all data we need    
    dic_data_covid = nx.get_node_attributes(G, name='data_covid')    
    # loop over time steps to get simu results and compute fit function
    for i in range(0,nbTimes):
        mytime = outputs[i][0]
        # not sure of what precision to take
        if (abs(mytime - time_date_target) <= 0.01) :
            print('[TARGET FOUND]: mytime, date_target = ' , mytime, time_date_target)
            x = outputs[i][1]
            total_confirmed = 0
            total_recovered = 0
            total_deaths    = 0
            for n in G.nodes():
                # Infected
                total_confirmed += dic_data_covid[n][params["dates_for_fit"][_count]]['Confirmed']
                total_recovered += dic_data_covid[n][params["dates_for_fit"][_count]]['Recovered']
                total_deaths    += dic_data_covid[n][params["dates_for_fit"][_count]]['Deaths']
            
            fitness[_count + 0] = (total_confirmed - x[1])
            fitness[_count + 1] = (total_recovered - x[2])
            fitness[_count + 2] = (total_deaths    - x[3])
            # get new date to target
            time_date_target = time_date_target + delta_days[_count]
            _count += 1
            # exist loop if we got all data we need
            if (_count == nbDates):
                break;

    print("fitness value = ", norm(fitness, ord=2) ** 2)
    return norm(fitness, ord=2) ** 2


def FitFuncSIRLDOneComp(paramsToFit, params, G): 
    print('FitFuncSIRLDOneComp')  
    # add parameters to fit into (that will be varied every call to the fitfunction)
    # to the dictionnary that contains all parameters
    params["k_SI"]  = float(paramsToFit["k_SI"])
    params["k_IR"]  = float(paramsToFit["k_IR"])
    params["k_IL"]  = float(paramsToFit["k_IL"])
    params["k_LD"]  = float(paramsToFit["k_LD"])
    params["L_ini"] = int(paramsToFit["L_ini"])
    params["I_ini"] = int(paramsToFit["I_ini"])
    params["alpha"] = float(paramsToFit["alpha"])
    
    print(paramsToFit)
    
    # transform to datetime format
    dates = [datetime.datetime.strptime(ts, "%Y-%m-%d") for ts in params["dates_for_fit"]] 
    time_start_date = datetime.datetime.strptime(params['time_start_day'], "%Y-%m-%d")        
    nbDates=len(params["dates_for_fit"])    
        
    # fill vector with relative increase from one date to the other    
    # dates[0] should NOT be the same than time_start_date otherwise 
    # delta_days[0] = 0 and the first fit data point is the intial values...
    # dates bmust be sorted in params['dates_for_fit']    
    delta_days = []
    delta_days.append((dates[0] - time_start_date).days) # transform to number of days
    for i in range(1,nbDates):
        delta_days.append((dates[i] - dates[i-1]).days) # transform to number of days
    
    # time_end must verify time_end = time_start + sum(delta_days)
    if params['time_end'] != params['time_start'] + sum(delta_days):        
        params['time_end'] = params['time_start'] + sum(delta_days)
        print("new time_end = " , params['time_end'])
    
    # init model and solve
    [params , G] = InitModels(params, G)
    outputs = SolveOneCompartment(params)
    nbTimes = len(outputs)
    fitness = np.zeros(3 * nbDates)
    time_date_target = params['time_start'] + delta_days[0]
    _count = 0
    # get all data we need    
    dic_data_covid = nx.get_node_attributes(G, name='data_covid')    
    # loop over time steps to get simu results and compute fit function
    x = np.zeros([nbTimes,params['nbStates']])
    _time = np.zeros(nbTimes)
    for i in range(0,nbTimes):        
        mytime = outputs[i][0]        
        _time[i] = mytime
        x[i,:] = outputs[i][1]
        # not sure of what precision to take
        if (abs(mytime - time_date_target) <= 0.01) :
            print('[TARGET FOUND]: mytime, date_target = ' , mytime, time_date_target)
            
            total_confirmed = 0
            total_recovered = 0
            total_deaths    = 0
            for n in G.nodes():
                # Infected       
                total_confirmed += dic_data_covid[n][params["dates_for_fit"][_count]]['Confirmed']
                total_recovered += dic_data_covid[n][params["dates_for_fit"][_count]]['Recovered']
                total_deaths    += dic_data_covid[n][params["dates_for_fit"][_count]]['Deaths']
            
            # compare confirmed to confirmed
            # do not add revoered to the fit function as it is not sure yet 
            # if there isn' t a delay like for the deaths
            # params["k_SI"] * params['alpha'] * x[1] 
            # tmp_funct = alpha * k_SI * S * I = beta I in marc lavielle
            tmp_func = params['alpha'] * params["k_SI"] * x[0:i,0] * x[0:i,1]            
            W_c = np.trapz(tmp_func, _time[0:i] ) # integral over time of I until mytime[i]
            print("W_c = ", W_c)
            fitness[_count + 0] = (total_confirmed - W_c)
            # fitness[_count + 1] = (total_recovered - x[2])
            fitness[_count + 2] = (total_deaths    - x[i,4])
            # get new date to target
            time_date_target = time_date_target + delta_days[_count]
            _count += 1
            # exist loop if we got all data we need
            if (_count == nbDates):
                break;
    
    params['num_iter_fit'] += 1
    
    print("num_iter_fit, fitness value = ", params['num_iter_fit'], norm(fitness, ord=2) ** 2)
    return norm(fitness, ord=2) ** 2

# parameters to fit AT THE MOMENT
# R_0 (between 2 and 4) --> this gives k_IS
# k_IR (between 1 / 10 and 1 / 20 days)
# k_ID (between (0.005 and 0.05)
def InitFit():
    print('InitFit')
    paramsToFit = lmfit.Parameters()
    # set parameters that need to be fitted
    paramsToFit.add("R_0", 3, min=1, max=5, vary=True)
    paramsToFit.add("k_IR", 0.0666, min=0.1, max=0.05, vary=True)
    paramsToFit.add("k_ID", 0.01, min=0.005, max=0.05, vary=True)
    return paramsToFit

def InitFitSIRLD():
    print('InitFitSIRLD')
    paramsToFit = lmfit.Parameters()    
    # set parameters that need to be fitted
    #beta between 0.01 and 1 
    paramsToFit.add("k_SI", value=1.4911195702987055e-9, min=1.4911195702987055e-10, max=1.4911195702987055e-8, vary=True)    
    #beta between 0.0001 and 0.001 to test quickly
    # paramsToFit.add("k_SI", value=1.4911195702987055e-9, min=1.4911195702987055e-12, max=1.4911195702987055e-11, vary=True)
    paramsToFit.add("k_IR", value=0.0666, min=0.1, max=0.05, vary=True)
    paramsToFit.add("k_IL", value=0.01, min=0.0001, max=1, vary=True)
    paramsToFit.add("k_LD", value=0.01, min=0.005, max=0.05, vary=True)
    paramsToFit.add("L_ini", value=400, min=100, max=1000, vary=True)
    paramsToFit.add("I_ini", value=400, min=5000, max=200000, vary=True)
    paramsToFit.add("alpha", value=0.1, min=0.01, max=1, vary=True)
    return paramsToFit


def LaunchFit(paramsToFit, params, G):
    print('LaunchFit')
    mini = lmfit.Minimizer(
        FitFuncSIRLDOneComp, paramsToFit, fcn_args=(params, G)
    )
    # mini = lmfit.Minimizer(
        # FitFuncOnGraph, paramsToFit, fcn_args=(params, G)
    # )
    results = mini.minimize(method="differential_evolution") #  plein de methods, least-square par default
    results.paramsToFit.pretty_print() #  nice display of fit results
    bic = results.bic #  example of info in results
    residuals = results.residuals #  example of info in results

    df = pd.DataFrame()  #  dataframe init
    df = df.append(results.paramsToFit.valuesdict(), ignore_index=True) #  tout le monde il est bien range
    
    params["R_0"]  = paramsToFit["R_0"]
    params["k_IR"] = paramsToFit["k_IR"]
    params["k_ID"] = paramsToFit["k_ID"]     
    
    return df , params
    
