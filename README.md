# ESSM_plotting
Author: Stanislaus J. Schymanski*

Luxembourg Institute of Science and Technology, Environmental Research and Innovation Department

*stanislaus.schymanski@list.lu

Copyright © 2023 Luxembourg Institute of Science and Technology. All Rights Reserved.

The code in this [repository](https://github.com/schymans/ESSM_plotting) is distributed under the terms of the [GNU General Public License 3.0](https://www.gnu.org/licenses/gpl-3.0.html) or any later version, unless stated otherwise. All other material is licensed under a [Creative Commons Attribution 4.0 International license](https://creativecommons.org/licenses/by/4.0/), unless stated otherwise.

## General purpose
This repository contains python scripts to plot symbolic expressions generated using the Python package [ESSM](https://essm.readthedocs.io). 

To use the functions in the file `plotting.py` within a Jupyter notebook, clone this repo as a submodule into the one where your jupyter notbook resides 
(e.g. `git submodule add https://github.com/schymans/ESSM_plotting.git ESSM_plotting`), 
then provide the relative path to the repo (e.g. `../ESSM_plotting/` and import the desired function:
```
import imp
path_plotting = '../ESSM_plotting/plotting.py`

# Importing plotting function
mod = imp.load_source('plotting', path_plotting)
plot_expr2 = getattr(mod, 'plot_expr2')
```

If you would like to import the functions e.g. into a project residing in [renkulab.io](https://renkulab.io/), execute in the parent folder of your renkulab project:
```
mkdir modules
renku dataset create --datadir modules/ESSM_plotting ESSM_plotting
renku --no-external-storage dataset add ESSM_plotting --source '*' git@github.com:schymans/ESSM_plotting.git

```
This will save the contents in `modules/ESSM_plotting`. 



