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
    
    print("params['k_SI']  = " , params["k_SI"] )
    print("params['k_IR'] = " , params["k_IR"] )
    print("params['k_IL'] = " , params["k_IL"] )
    print("params['k_LD'] = " , params["k_LD"] )
    print("params['L_ini'] = " , params["L_ini"])
    print("params['I_ini'] = " , params["I_ini"])
    print("params['alpha'] = " , params["alpha"])    
        
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
        W_c = 0.
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
            
            tmp_func = params['alpha'] * params["k_SI"] * x[0:i+1,0] * x[0:i+1,1]
            W_c = np.trapz(tmp_func, _time[0:i+1] ) # integral over time of I until _time[i]  
            # add integral from 0 to t_0 from data so we can compare
            W_c += params['total_covid_data'][params['time_start_day']]['Confirmed']            
            fitness[_count + 0] = (total_confirmed - W_c)
            # recovered in model is a fraction of recovered measured
            fitness[_count + 1] = (total_recovered - params['alpha'] * x[i,2])
            fitness[_count + 2] = (total_deaths    - x[i,4])
            # get new date to target
            time_date_target = time_date_target + delta_days[_count]
            _count += 1
            # exist loop if we got all data we need
            if (_count == nbDates):
                break;
    
    params['num_iter_fit'] += 1    
    fitness_L2_square = norm(fitness, ord=2) ** 2
    print("num_iter_fit, fitness value = ", params['num_iter_fit'], fitness_L2_square)

    return fitness_L2_square

def FitFuncSEIRLDOneComp(paramsToFit, params, G): 
    print('FitFuncSEIRLDOneComp')  
    # add parameters to fit into (that will be varied every call to the fitfunction)
    # to the dictionnary that contains all parameters
    params["k_SE"]  = float(paramsToFit["k_SE"])
    params["k_EI"]  = float(paramsToFit["k_EI"])
    params["k_IR"]  = float(paramsToFit["k_IR"])
    params["k_IL"]  = float(paramsToFit["k_IL"])
    params["k_LD"]  = float(paramsToFit["k_LD"])
    params["k_LR"]  = float(paramsToFit["k_LR"])
    params["alpha"] = float(paramsToFit["alpha"])
    params["E_ini"] = int(paramsToFit["E_ini"])
    params["I_ini"] = int(paramsToFit["I_ini"])
    params["L_ini"] = int(paramsToFit["L_ini"])
    
    print("params['k_SE']  = " , params['k_SE'] )
    print("params['k_EI']  = " , params['k_EI'] )
    print("params['k_IR']  = " , params['k_IR'] )
    print("params['k_IL']  = " , params['k_IL'] )
    print("params['k_LD']  = " , params['k_LD'] )
    print("params['k_LR']  = " , params['k_LR'] )
    print("params['alpha'] = " , params['alpha'])    
    print("params['E_ini'] = " , params['E_ini'])    
    print("params['I_ini'] = " , params['I_ini'])    
    print("params['L_ini'] = " , params['L_ini'])    
        
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
    # get all data we need    
    dic_data_covid = nx.get_node_attributes(G, name='data_covid')    
    dic_total_pop  = nx.get_node_attributes(G, name='total_pop')
    
    # get total population
    N_pop = 0 
    for n in G.nodes():
        N_pop += dic_total_pop[n] 
    
    # loop over time steps to get simu results and compute fit function
    x = np.zeros([nbTimes,params['nbStates']])
    _time = np.zeros(nbTimes)
    _count = 0
    for i in range(0,nbTimes):        
        mytime = outputs[i][0]        
        _time[i] = mytime
        x[i,:] = outputs[i][1]        
        # not sure of what precision to take
        W_c = 0.
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
            
            tmp_func = params['alpha'] * params["k_SE"] * x[0:i+1,0] * x[0:i+1,2] / N_pop
            # print("x[0:i+1,0, x[0:i+1,2], N_pop, tmp_func " , x[0:i+1,0] , x[0:i+1,2], N_pop, tmp_func)
            W_c = np.trapz(tmp_func, _time[0:i+1] ) # integral over time of I until _time[i]  
            # add integral from 0 to t_0 from data so we can compare
            W_c += params['total_covid_data'][params['time_start_day']]['Confirmed']
            fitness[_count + 0] = (total_confirmed - W_c)
            # recovered in model is a fraction of recovered measured
            fitness[_count + 1] = (total_recovered - params['alpha'] * x[i,3])
            fitness[_count + 2] = (total_deaths    - x[i,5])
            # print(total_confirmed, total_recovered, total_deaths)
            # print(x[i,:], W_c)
            # print("\n")
            
            # get new date to target
            time_date_target = time_date_target + delta_days[_count]
            _count += 1
            # exist loop if we got all data we need
            if (_count == nbDates):
                break;
            
    if not(_count == nbDates):
        print("[ERROR] : we didn't find all dates in model to compute fitness function".)
    # print('fitness = ' , fitness)
    params['num_iter_fit'] += 1    
    fitness_L2_square = norm(fitness, ord=2) ** 2
    print("num_iter_fit, fitness value = ", params['num_iter_fit'], fitness_L2_square)

    return fitness_L2_square

def InitFitSIRLD():
    print('InitFitSIRLD')
    paramsToFit = lmfit.Parameters()    
    # set parameters that need to be fitted
    #beta between 0.01 and 1 
    paramsToFit.add("k_SI", value=1.5e-10, min=1e-10, max=1e-9, vary=True)    
    #beta between 0.0001 and 0.001 to test quickly
    # paramsToFit.add("k_SI", value=1.4911195702987055e-9, min=1.4911195702987055e-12, max=1.4911195702987055e-11, vary=True)
    paramsToFit.add("k_IR", value=0.0666, min=0.0001, max=0.1, vary=True)
    paramsToFit.add("k_IL", value=0.01, min=0.0001, max=0.1, vary=True)
    paramsToFit.add("k_LD", value=0.01, min=0.0001, max=0.1, vary=True)
    paramsToFit.add("L_ini", value=400, min=100, max=1000, vary=True)
    paramsToFit.add("I_ini", value=400, min=5000, max=200000, vary=True)
    paramsToFit.add("alpha", value=0.1, min=0.01, max=1, vary=True)
    return paramsToFit

def InitFitSEIRLD():
    print('InitFitSEIRLD')
    paramsToFit = lmfit.Parameters()    
    # set parameters that need to be fitted
    # kinetic parameters  
    # we actually set k_SE_tilde = k_SE * N. in InitModels k_SE is divided by N
    paramsToFit.add("k_SE", value=0.67063703, min=0.1, max=1, vary=True)
    paramsToFit.add("k_EI", value=0.0666, min=0.0001, max=0.1, vary=True)
    paramsToFit.add("k_IR", value=0.01, min=0.0001, max=0.1, vary=True)
    paramsToFit.add("k_IL", value=0.01, min=0.0001, max=0.1, vary=True)
    paramsToFit.add("k_LD", value=400, min=100, max=1000, vary=True)
    paramsToFit.add("k_LR", value=400, min=5000, max=200000, vary=True)
    # fraction of infected in population I_c = alpha * I 
    paramsToFit.add("alpha", value=0.1, min=0.01, max=1, vary=True)
    # initial values
    paramsToFit.add("E_ini", value=400, min=1, max=4000, vary=True)
    paramsToFit.add("I_ini", value=400, min=1, max=2000, vary=True)
    paramsToFit.add("L_ini", value=400, min=1, max=1000, vary=True)
    
    return paramsToFit

def InitFit(params):
    print('InitFit')
    paramsToFit = lmfit.Parameters()    
    if (params["model_name"] == "RealSIRLDWellMixed"):
        paramsToFit = InitFitSIRLD()
    elif (params["model_name"] == "RealSEIRLDWellMixed"):
        paramsToFit = InitFitSEIRLD()
    else: 
        print('[ERROR] : unknown model name: ',params["model_name"] ,' in InitFit function model ')
    
    return paramsToFit


def LaunchFit(paramsToFit, params, G):
    print('LaunchFit')    
    if (params["model_name"] == "RealSIRLDWellMixed"):
        mini = lmfit.Minimizer(
        FitFuncSIRLDOneComp, paramsToFit, fcn_args=(params, G))
    elif (params["model_name"] == "RealSEIRLDWellMixed"):
        mini = lmfit.Minimizer(
        FitFuncSEIRLDOneComp, paramsToFit, fcn_args=(params, G))    
    else: 
        print('[ERROR] : unknown model name: ',params["model_name"] ,' in LaunchFit function model ')
    
    results = mini.minimize(method=str(params['method'])) #  plein de methods, least-square par default
    results.paramsToFit.pretty_print() #  nice display of fit results
    bic = results.bic #  example of info in results
    residuals = results.residuals #  example of info in results

    df = pd.DataFrame()  #  dataframe init
    df = df.append(results.paramsToFit.valuesdict(), ignore_index=True) #  tout le monde il est bien range
    
    return df , params
    
