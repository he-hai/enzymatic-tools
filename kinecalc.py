#%%
import numpy as np 
from lmfit import minimize, Parameters
import pandas as pd 

def MMfunc(params, E, S, v=None):
    '''Michaelis–Menten kinetics, first-order  
    :params kcat: turnover number
    :params Km: Michaelis constant
    :params S: substrate concentration
    :params v: reaction rate
    :params E: enzyme concentration
    '''
    kcat = params['kcat']
    Km = params['Km']

    model = kcat * E * S / (Km + S)
    if v is None:
        return model
    return model - v

params = Parameters()
params.add('kcat', value=10, min=0)  # in s-1
params.add('Km', value=1e-4, min=0)  # in M


#%%
S = np.array([0.5, 0.25, 0.125, 0.05, 5]) * 1e-3
v = np.array([0.249, 0.126, 0.058, 0.027, 0.357]) * 1e-3 / 60
E = 1e-6

kinetics = minimize(MMfunc, params, args=(E, S,), kws={'v': v})
print('kcat: %.2f s-1' % kinetics.params['kcat'].value)
# print(kinetics.params['vmax'].value)
print('Km: %.2e M' % kinetics.params['Km'].value)
print(kinetics.chisqr)
# %%
