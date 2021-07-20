#%%
import numpy as np 
from lmfit import minimize, Parameters
# import pandas as pd 
import pint 
import matplotlib.pyplot as plt 

ureg = pint.UnitRegistry()
Q_ = ureg.Quantity
ureg.default_format = "~P"
ureg.context('chemistry')

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

def kinetics_calc(E, S, v):
    kinetics_fit = minimize(MMfunc, params, args=(E, S,), kws={'v': v})
    Km = kinetics_fit.params['Km'].value
    kcat = kinetics_fit.params['kcat'].value
    chisqr = kinetics_fit.chisqr
    return Km, kcat, chisqr

def kinetics_report(
    E, S, v,
    enzyme=None, # name
    substrate=None,
    **kwg
):
    Km, kcat, chisqr = kinetics_calc(E, S, v)

    fig, ax = plt.subplots()
    x_up = 1.2 * S.max()
    x = np.arange(0,x_up,0.01*x_up)
    y = kcat * E * x / (Km + x)
    ax.plot(x, y,c='blue')
    ax.scatter(S, v, c='orange')

    ax.hlines(kcat*E, 0, x_up, colors='grey', linestyles='--')
    ax.vlines(Km, 0, kcat*E/2,colors='grey',linestyles='--')
    ax.text(1.2*Km, kcat*E/4,
       r'$K_m$ ' f'{Km:.2e}\n' 
       r'$k_{cat}$ ' f'{kcat:.2f}\n' 
       r'$\chi^2$ ' f'{chisqr:.3e}',
       )
    ax.title.set_text(f'{enzyme} - {substrate}')
    ax.set(xlabel=f'[{substrate}]', ylabel='activity')
    ax.set_xlim(xmin=0)
    ax.set_ylim(ymin=0)
    
    fig.savefig(f'{enzyme}-{substrate} kinetics', dpi=300)
    fig.show()

#%%
Enzyme = 'enzyme name'
Substrate = 'Substrate'
S = np.array([0.5, 0.25, 0.125, 0.05, 5]) * 1e-3
v = np.array([0.249, 0.126, 0.058, 0.027, 0.357]) * 1e-3 / 60
E = 1e-6

kinetics_report(E, S, v, Enzyme, Substrate)