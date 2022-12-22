###****
# For calculating kinetics with substrate inhibition
###****
#%%
import numpy as np 
from lmfit import Model
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
# plt.rcParams.update({'mathtext.default': 'regular'})
plt.rcParams['axes.axisbelow'] = True
plt.rcParams['savefig.bbox']='tight'

'''Main function/calculations'''
def MMfunc(S, kcat, Km, Ki):
    '''Substrate inibition model 
    :params S: substrate concentration
    :params kcat: turnover number
    :params Km: Michaelis-Menten constant
    :params Ki: dissociation constant
    '''
    return kcat * S / (Km + S * (1 + S / Ki))

model = Model(MMfunc, nan_policy='omit')
model.set_param_hint('kcat',value=10,min=0)  # in s-1
model.set_param_hint('Km',value=1e-4,min=0)  # in M
model.set_param_hint('Ki',value=1e-2,min=0)  # in M

def kinetics_calc(S, t):
    result = model.fit(t, S=S)
    print(result.fit_report())
    Km = Q_(
        result.params['Km'].value,'M'
    ).plus_minus(
        result.params['Km'].stderr
    )
    kcat = Q_(
        result.params['kcat'].value,'1/s'
    ).plus_minus(
        result.params['kcat'].stderr
    )
    Ki = Q_(
        result.params['Ki'].value,'M'
    ).plus_minus(
        result.params['Ki'].stderr
    )
    chisqr = result.chisqr
    return Km, kcat, Ki, chisqr

def kinetics_report(
    S, t,
    enzyme=None, # name
    substrate=None,
    lb = Q_(0,'M'),
    ub = Q_(np.Infinity,'M'), 
    **kwg
):
    # E_ = E.to('M').magnitude
    indx = np.where((S>=lb) & (S<=ub))
    S_ = S[indx].to('M').magnitude
    S_u = S.units
    # v_ = v.to('M/s').magnitude
    t_ = t[indx].magnitude
    
    Km, kcat, Ki, chisqr = kinetics_calc(S_, t_)
    cateff = kcat/Km
    Km = Km.to(S_u)
    Ki = Ki.to(S_u)

    fig, ax = plt.subplots(figsize=(8,6))
    x_up = 1.2 * S.max()
    x = np.linspace(0,x_up,100)
    y = MMfunc(x, kcat.value, Km.value, Ki.value)
    ax.scatter(S, t, c='orange')
    ax.plot(x, y, c='blue')

    ax.hlines(kcat.value, 0, x_up, colors='grey', linestyles='--')
    ax.vlines(Km.value, 0, kcat.value/2,colors='grey',linestyles='--')
    ax.text(2*Km.value, kcat.value/5,
       '$\chi^2$: ' f'{chisqr:.3e}\n'
       '$K_m$: ' f'{Km.magnitude:.2f} {Km.units:~P}\n'
       '$k_{cat}$: ' f'{kcat.magnitude:.2f} {kcat.units:~P}\n'
       '$K_i$: ' f'{Ki.magnitude:.2f} {Ki.units:~P}\n'
       '$k_{cat} / K_m$: ' f'{cateff.magnitude:.2e} {cateff.units:~P}' ,
       fontsize=15,
       )
    ax.title.set_text(f'{enzyme} - {substrate}')
    ax.set(xlabel=f'[{substrate}] ({S_u})', 
           ylabel='$v_0/[E]\ (s^{-1}$)')
    ax.set_xlim(xmin=0)
    ax.set_ylim(ymin=0)
    
    fig.savefig(f'{enzyme}-{substrate} kinetics', dpi=300)
    fig.show()

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

S = Q_(df['substrate'].to_numpy(),'mM')
S = S.to('uM')  # change to unit for plot if different
# ub = Q_(0.1,'mM')  # upper bound for model calculations, included
# lb = Q_(0.005, 'mM')  # lower bound for model calculations, included

v = Q_(df['slope'].to_numpy(),'1/min')/(Enadph*l)  # initial velocity
t = (v/E).to('1/s')  # turn over

# check with Lineaweaver Burk Plot 
sns.regplot(
    x=1/S,y=1/v,
)

kinetics_report(S, t, Enzyme, Substrate)
