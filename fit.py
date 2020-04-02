import lmfit
import numpy as np
from numpy.linalg import norm
import pandas as pd

#  def function(params, dependant variables):
def fitFunc(params, t, Cpv, Clv, options):

    initmodel
    [outputs , G] = solve 
    get closest time 
    nbDates=len(options["dates_for_fit"])
    fitness = np .zeros(nbDates)
    
    nbTimes = len(outputs)
    dates_int=[int(d.replace('-','')) for d in options["dates_for_fit"]]
    dates_int.sort()
    time_start = options["time_start"]
    delta_dates = [dates_int[i] - dates_int[i-1] for i in range(1,len(dates_int))]
    time_target = dates_int[0]
    for i in range(0,nbTimes):
        mytime = outputs[i][0]
        
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

