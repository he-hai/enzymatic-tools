This repository contains scripts for enzyme kinetic parameters calculation.  

`kinecalc.py` is for Michaelis–Menten kinetics $v=\frac{Vmax*[S]}{Km+[S]}$.  

And `kinecalc_subI.py` is for kinetics with substrate inhibition $v=\frac{Vmax*[S]}{Km+[S]+\frac{[S]^2}{Ki}}$.  

The scripts use [lmfit](https://lmfit.github.io/lmfit-py/model.html) for non-linear curve fitting, [pint](https://pint.readthedocs.io/en/stable/) for units handling. 

The scripts have been tested in an environment with following packages:  
- python 3.8.2
- numpy 1.23.4
- pandas 1.3.4
- matplotlib 3.5.1
- seaborn 0.11.2
- lmfit 1.0.0
- openpyxl 3.0.6
- scipy 1.4.1
- pint 0.17

### Usage
The scripts take directly raw data, such as slope (in the unit of _1/min_), enzyme concentration (in molar concentrations, _$\mu$M_ ), extinction coefficient. To use the script, please adapt these input directly in the script. Raw data can be in the format as in the [Demo-data.xlsx](./Demo-data.xlsx). `pint` will handle all the unit conversions, so please adapt the units if necessary. 