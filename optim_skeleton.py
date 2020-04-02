import lmfit
import numpy as np
from numpy.linalg import norm
import pandas as pd

#  def function(params, dependant variables):
def pk_1cCP(params, t, Cpv, Clv):
    k_h = params["k_h"]
    F = params["F"]

    inInput = F * Cpv

    C = np.convolve(np.exp(-k_h * t), inInput * 1)[: len(t)]

    return norm(C - Clv, ord=2) ** 2


params = lmfit.Parameters()
params.add("F", 0.01, min=0, max=15, vary=False)
params.add("k_h", 1, min=0, max=0.5, vary=True)


mini = lmfit.Minimizer(
    pk_1cCP, params, fcn_args=(t_clv, Cpvi, Clvi)
)
results = mini.minimize(method="differential_evolution") #  plein de methods, least-square par default
results.params.pretty_print() #  nice display of fit results
bic = results.bic #  example of info in results
residuals = results.residuals #  example of info in results

df = pd.DataFrame()  #  dataframe init
df = df.append(results.params.valuesdict(), ignore_index=True) #  tout le monde il est bien range

