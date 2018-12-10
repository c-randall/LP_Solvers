"""
Created on Mon Oct 29 15:27:08 2018

Author: 
    Corey R. Randall

Summary:
    Phase I of Simplex Method coded in Python for MEGN 586. Work completed at 
    Colorado School of Mines during the Fall 2018 semester.

"""

""" Import Useful Modules """
"-----------------------------------------------------------------------------"
import numpy as np

""" User Inputs """
"-----------------------------------------------------------------------------"
# Objective Function:
objective = 'maximize'
c_coeff = [13, 16, 16, 14, 39] 

# Constraints:
constraint = {}
constraint[1] = ['time_0', [11, 53, 5, 5, 29], '<=', 40]
constraint[2] = ['time_1', [3, 6, 5, 1, 34], '<=', 20]
for i in range(5):
    frac_coeff = np.zeros(5)
    frac_coeff[i] = 1
    constraint[i+3] = ['frac_'+str(i+1), frac_coeff.tolist(), '<=', 1]

# Toggles/options for sensitivity and incoming basis selections:
sensitivity = 'on'   # runs a RHS sensitivity analysis ('on' or 'off')
incoming = 'last'    # incoming varialbe ('first' or 'last') neg c_bar value

# Tolerance Conditions: 
# Makes the (abs(x) <= tolerance) equal 0 in c_bar calculations
tolerance = 1e-10
decimals = 3          # up to this many decimals will be printed in output

""" Report User Input Errors """
"-----------------------------------------------------------------------------"
# Fill this in with errors to check lengths of variables... constraints... etc.

""" Create Pointers for Later Modifications """
"-----------------------------------------------------------------------------"
pointer = {}
pointer['name'] = 0
pointer['coeff'] = 1
pointer['inequality'] = 2
pointer['b-value'] = 3

""" Run Standard Form Conversion """
"-----------------------------------------------------------------------------"
# Build dictionary to pass information to subfunctions:
user_inputs = {}
user_inputs['objective'] = objective
user_inputs['c_coeff'] = c_coeff
user_inputs['tolerance'] = tolerance
user_inputs['pointer'] = pointer
user_inputs['incoming'] = incoming
user_inputs['dec'] = decimals

from Standard_Form import standard_form
conversion = standard_form(user_inputs, constraint)

""" Run Phase I """
"-----------------------------------------------------------------------------"
from Phase_1 import phase_1
initial_solution = phase_1(user_inputs, conversion)

""" Run Phase II """
"-----------------------------------------------------------------------------"
from Phase_2 import phase_2
solution = phase_2(user_inputs, conversion, initial_solution)

""" If specified, perform a RHS Sensitivity Analysis """
"-----------------------------------------------------------------------------"
if all([sensitivity == 'on', initial_solution['feasibility'] == 'feasible',
        solution['bounded'] == 'yes']):
    from RHS_Sensitivity import rhs_sensitivity
    rhs_sensitivity_report = rhs_sensitivity(user_inputs, constraint, 
                                             conversion, solution)
