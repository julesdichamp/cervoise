import lmfit
import numpy as np
from numpy.linalg import norm
import pandas as pd
from models import InitModels, SolveOnGraph
from datetime import date

#  def function(params, dependant variables):
def fitFunc(params, t, otherParams, G):

    [otherParams , G] = InitModels(otherParams, G)
    [outputs , G] = SolveOnGraph(otherParams, G)
    nbDates=len(options["dates_for_fit"])
    fitness = np.zeros(nbDates)
    
    nbTimes = len(outputs)    
    dates = [datetime.datetime.strptime(ts, "%Y-%m-%d") for ts in options["dates_for_fit"]]
    
    time_start_date = datetime.datetime.striptime(options['time_start_date'], "%Y-%m-%d")
    delta_days = []
    delta_days[0] = dates[0] - time_start_date
    for i in range(1,nbDates):
        delta_days[i] = dates[i] - dates[i-1]
        
    date_target = time_start_date + delta_days[0]
    _count = 0
    for i in range(0,nbTimes):
        mytime = outputs[i][0]
        # not sure of what precision to take
        if (abs(mytime - date_target)) <= 0.001:
            date_target = date_target + delta_days[_count]
            _count += 1
        dic_x = outputs[i][1]    
        
    C = np.convolve(np.exp(-k_h * t), inInput * 1)[: len(t)]
    

    return norm(C - Clv, ord=2) ** 2


params = lmfit.Parameters()
# set parameters that need to be fitted
# params.add("F", 0.01, min=0, max=15, vary=False)
# params.add("k_h", 1, min=0, max=0.5, vary=True)


mini = lmfit.Minimizer(
    fitFunc, params, fcn_args=(t_clv, Cpvi, Clvi)
)
results = mini.minimize(method="differential_evolution") #  plein de methods, least-square par default
results.params.pretty_print() #  nice display of fit results
bic = results.bic #  example of info in results
residuals = results.residuals #  example of info in results

df = pd.DataFrame()  #  dataframe init
df = df.append(results.params.valuesdict(), ignore_index=True) #  tout le monde il est bien range

