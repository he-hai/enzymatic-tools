#%%
import numpy as np 
from lmfit import minimize, Parameters
import pandas as pd 
import pint 
import matplotlib.pyplot as plt 
import seaborn as sns 

ureg = pint.UnitRegistry()
Q_ = ureg.Quantity
ureg.default_format = "~P"
ureg.context('chemistry')
ureg.setup_matplotlib(True)

plt.rcParams["font.family"] = "Arial"
plt.rcParams['ytick.labelsize'] = 15
plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['axes.titlesize'] = 24
plt.rcParams['axes.labelsize'] = 20
plt.rcParams['legend.title_fontsize'] = 15
plt.rcParams['legend.fontsize'] = 13
plt.rcParams.update({'mathtext.default': 'regular'})
plt.rcParams['axes.axisbelow'] = True
plt.rcParams['savefig.bbox']='tight'
#%%
'''Region for inputs'''
Enzyme = 'enzyme name'
Substrate = 'Substrate'
# S = Q_(np.array([0.5, 0.25, 0.125, 0.05, 5]),'mM')
# v = Q_(np.array([0.249, 0.126, 0.058, 0.027, 0.357]),'mM/min')
E = Q_(0.02,'uM')

Enadph = Q_(6220,'1/(M*cm)')  # M-1 cm-1
l = Q_(1,'cm')  # cm

df = pd.read_excel(
    'Demo-data.xlsx',sheet_name='Sheet1',
    engine='openpyxl', index_col=0,
)

df= df.stack().reset_index()
df.rename(columns={'level_1':'replicates',0:'slope'},inplace=True)

S = Q_(df['substrate'].to_numpy(),'mM').to('uM')
v = Q_(df['slope'].to_numpy(),'1/min')/(Enadph*l)

# check with Lineaweaver Burk Plot 
sns.regplot(
    x=1/S,y=1/v,
)
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
    kinetics_fit = minimize(
        MMfunc, params, args=(E, S,), kws={'v': v},
        nan_policy='omit'
    )
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

    fig, ax = plt.subplots(figsize=(8,6))
    x_up = 1.2 * S.max()
    x = np.linspace(0,x_up,100)
    y = (kcat * E * x / (Km + x)).to(ureg('uM/min'))
    ax.scatter(S, v, c='orange')
    ax.plot(x, y,c='blue')

    ax.hlines(kcat*E, 0, x_up, colors='grey', linestyles='--')
    ax.vlines(Km, 0, kcat*E/2,colors='grey',linestyles='--')
    ax.text(1.2*Km, kcat*E/4,
       '$\chi^2$: ' f'{chisqr:.3e}\n' 
       '$K_m$: ' f'{Km.to(S_u):.2f~P}\n' 
       '$k_{cat}$: ' f'{kcat:.2f~P}\n' 
       '$k_{cat} / K_m$: ' f'{kcat/Km:.2e~P}',
       )
    ax.title.set_text(f'{enzyme} - {substrate}')
    ax.set(xlabel=f'[{substrate}] ({S_u})', 
           ylabel='Activity ($\mu$M/min)')
    ax.set_xlim(xmin=0)
    ax.set_ylim(ymin=0)
    
    fig.savefig(f'{enzyme}-{substrate} kinetics', dpi=300)
    fig.show()

kinetics_report(E, S, v, Enzyme, Substrate)