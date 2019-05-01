"""
Created on Mon Oct 29 15:27:08 2018

Author: 
    Corey R. Randall

Summary:
    Input and runner file for Primal Simplex Method. Inputs can also be loaded
    from a separate file for larger problems. An example of a script for a 
    larger problem is shown in pv_batt_sizing.py and pv_batt_sizing_sub.py.

"""

""" Import Useful Modules """
"-----------------------------------------------------------------------------"
import numpy as np
import os; cwd = os.getcwd();

""" User Inputs """
"-----------------------------------------------------------------------------"
# Objective Function: --> ignore due to script import below
objective = ''
c_coeff = [] 

# Constraints: --> ignore due to script import below
constraint = {}
constraint[1] = []

" Alternative to filling in script obj. func. and constraints --> import "
os.chdir(cwd + '/Problem_Scripts')
from pv_batt_sizing_sub import v, objective, c_coeff, constraint

# Toggles/options:
sensitivity = 'off'   # runs a RHS sensitivity analysis ('on' or 'off')
pricing = 'most'      # pricing scheme for phase 2 ('most' or 'bland')
incoming = 'first'    # incoming varialbe ('first' or 'last') for bland's
problem = 'primal'    # problem that method operates on ('primal' or 'dual')
method = 'primal'     # algorithm method ('primal' or 'dual') Simplex

# Tolerance Conditions: 
tolerance = 1e-6      # treat any abs(c_bar) <= tolerance as 0
decimals = 3          # maximum number of decimals printed in output

# Display Variables:
var_names = ['x_B', 'x_PV']     

""" End of user inputs - do not edit anything below this line """
"-----------------------------------------------------------------------------"
###############################################################################
###############################################################################
###############################################################################

""" Execute LP_Simplex_Runner.py with user inputs from above  """
"-----------------------------------------------------------------------------"
if __name__ == '__main__':
    exec(open(cwd + '/Simplex_Files/LP_Simplex_Runner.py').read())
    os.chdir(cwd)
