import lmfit
import numpy as np
from numpy.linalg import norm
import pandas as pd
from models import InitModels, SolveOnGraph
# from datetime import date
import datetime
import networkx as nx

def FitFunc(paramsToFit, params, G): 
    print('FitFunc')  
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
    
    # time_end must verify time_end = time_start + delta_days[nbDates -1]
    if params['time_end'] != params['time_start'] + delta_days[nbDates-1]:
        params['time_end'] = params['time_end'] + delta_days[nbDates-1]
        print(params['time_end'])
    
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


def LaunchFit(paramsToFit, params, G):
    print('LaunchFit')
    mini = lmfit.Minimizer(
        FitFunc, paramsToFit, fcn_args=(params, G)
    )
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
    
