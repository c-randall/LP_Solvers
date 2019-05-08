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
# objective = 'minimize' or 'maximize'
# c_coeff = [coeff1, coeff2, ...] 

# Constraints: --> ignore due to script import below
# constraint = {}
# constraint[1] = ['name1', [coeff1, coeff2, ...], '<=' '>=' or '=', b]
# constraint[2] = [...]

# Alternative to in-file inputs for obj. func. and constraints --> import
" Import a problem script from a file in the same directory as LP_User_Inputs "
# from file_name import v, objective, c_coeff, constraint

" Import a problem script from a folder named Problem_Scripts "
# os.chdir(cwd + '/Problem_Scripts')
# from file_name import v, objective, c_coeff, constraint
# os.chdir(cwd)

# Settings/options:
sensitivity = 'off'   # runs a RHS sensitivity analysis ('on' or 'off')
pricing = 'most'      # pricing scheme for phase 2 ('most' or 'bland')
incoming = 'first'    # incoming varialbe ('first' or 'last') for bland's
problem = 'primal'    # problem that method operates on ('primal' or 'dual')
method = 'dual'       # algorithm method ('primal' or 'dual') Simplex
verbose = 'off'       # print each iteration's objective value ('on' or 'off')

# Tolerance Conditions: 
tolerance = 1e-6      # abs(data or calcs) <= tolerance as 0 (use scientific)
decimals = 3          # maximum number of decimals printed in output

# Display Variables:
var_names = None  

""" End of user inputs - do not edit anything below this line """
"-----------------------------------------------------------------------------"
###############################################################################
if __name__ == '__main__':
    exec(open(cwd + '/Simplex_Files/LP_Simplex_Runner.py').read())
