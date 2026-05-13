# ESSM_plotting
Author: Stanislaus J. Schymanski*

Luxembourg Institute of Science and Technology, Environmental Research and Innovation Department

*stanislaus.schymanski@list.lu

Copyright © 2023 Luxembourg Institute of Science and Technology. All Rights Reserved.

The code in this [repository](https://github.com/schymans/ESSM_plotting) is distributed under the terms of the [GNU General Public License 3.0](https://www.gnu.org/licenses/gpl-3.0.html) or any later version, unless stated otherwise. All other material is licensed under a [Creative Commons Attribution 4.0 International license](https://creativecommons.org/licenses/by/4.0/), unless stated otherwise.

## General purpose
This repository contains python scripts to plot symbolic expressions generated using the Python package [ESSM](https://essm.readthedocs.io). 

To use the functions in the file `plotting.py` within a Jupyter notebook, clone this repo as a submodule into the one where your jupyter notbook resides, e.g.:
```
mkdir submodules
cd submodules
git submodule add https://github.com/schymans/ESSM_plotting.git ESSM_plotting
```
Then provide the relative path to the repo (e.g. `../ESSM_plotting/` and import the desired function:
```
import imp
path_plotting = '../ESSM_plotting/plotting.py`

# Importing plotting function
mod = imp.load_source('plotting', path_plotting)
plot_expr2 = getattr(mod, 'plot_expr2')
```
If you or someone else creates a fresh clone
of your repo, the submodules can be pulled in with the following command run in the base folder:
```
git submodule init
git submodule update
```
At any stage, you can update all submodules in your repo with a single command:
```
git submodule update --remote --merge
```




