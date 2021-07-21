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
ureg.setup_matplotlib(True)

#%%
'''Region for inputs'''
Enzyme = 'enzyme name'
Substrate = 'Substrate'
S = Q_(np.array([0.5, 0.25, 0.125, 0.05, 5]),'mM')
v = Q_(np.array([0.249, 0.126, 0.058, 0.027, 0.357]),'mM/min')
E = Q_(1,'uM')

#%%
'''Main function/calculations'''
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
    Km = Q_(kinetics_fit.params['Km'].value,'M')
    kcat = Q_(kinetics_fit.params['kcat'].value,'1/s')
    chisqr = kinetics_fit.chisqr
    return Km, kcat, chisqr

def kinetics_report(
    E, S, v,
    enzyme=None, # name
    substrate=None,
    **kwg
):
    E_ = E.to('M').magnitude
    S_ = S.to('M').magnitude
    S_u = S.units
    v_ = v.to('M/s').magnitude
    
    Km, kcat, chisqr = kinetics_calc(E_, S_, v_)

    fig, ax = plt.subplots()
    x_up = 1.2 * S.max()
    x = np.linspace(0,x_up,100)
    y = (kcat * E * x / (Km + x)).to(ureg('uM/min'))
    ax.plot(x, y,c='blue')
    ax.scatter(S, v, c='orange')

    ax.hlines(kcat*E, 0, x_up, colors='grey', linestyles='--')
    ax.vlines(Km, 0, kcat*E/2,colors='grey',linestyles='--')
    ax.text(1.2*Km, kcat*E/4,
       r'$\chi^2$: ' f'{chisqr:.3e}\n' 
       r'$K_m$: ' f'{Km.to(S_u):.2f~P}\n' 
       r'$k_{cat}$: ' f'{kcat:.2f~P}\n' 
       r'$k_{cat} / K_m$: ' f'{kcat/Km:.2e~P}',
       )
    ax.title.set_text(f'{enzyme} - {substrate}')
    ax.set(xlabel=f'[{substrate}] ({S_u})', 
           ylabel=r'Activity ($\mu$M/min)')
    ax.set_xlim(xmin=0)
    ax.set_ylim(ymin=0)
    
    fig.savefig(f'{enzyme}-{substrate} kinetics', dpi=300)
    fig.show()

kinetics_report(E, S, v, Enzyme, Substrate)