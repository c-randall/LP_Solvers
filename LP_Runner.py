"""
Created on Mon Oct 29 15:27:08 2018

Author: 
    Corey R. Randall

Summary:
    Phase I of Simplex Method coded in Python for MEGN 586. Work completed at 
    Colorado School of Mines during the Fall 2018 semester.

"""

""" User Inputs """
"-----------------------------------------------------------------------------"
# Objective Function: --> ignore due to script import below
 

# Constraints: --> ignore due to script import below


" Alternative to filling in script obj. func. and constraints --> import "
from pv_batt_sizing import v, objective, c_coeff, constraint

# Toggles/options for sensitivity and incoming basis selections:
sensitivity = 'off'   # runs a RHS sensitivity analysis ('on' or 'off')
pricing = 'steep'     # pricing scheme for phase 2 ('steep' or 'bland')
incoming = 'first'    # incoming varialbe ('first' or 'last') for bland's

# Tolerance Conditions: 
tolerance = 1e-6      # treat any abs(c_bar) <= tolerance as 0
decimals = 3          # up to this many decimals will be printed in output

# Display Variables:
var_names = ['x_B', 'x_PV']     # list of variables to display after solution

""" Create Pointers for Later Modifications """
"-----------------------------------------------------------------------------"
pointer = {}
pointer['name'] = 0
pointer['coeff'] = 1
pointer['inequality'] = 2
pointer['b-value'] = 3

""" Run Standard From Conversion """
"-----------------------------------------------------------------------------"
# Build dictionary to pass information to subfunctions:
user_inputs = {}
user_inputs['objective'] = objective
user_inputs['c_coeff'] = c_coeff
user_inputs['tolerance'] = tolerance
user_inputs['pointer'] = pointer
user_inputs['pricing'] = pricing
user_inputs['incoming'] = incoming
user_inputs['dec'] = decimals

from Standard_Form import standard_form
conversion = standard_form(user_inputs, constraint)

""" Pre-Processing """
"-----------------------------------------------------------------------------"
# Fill this in with subfunction that simplifies and checks for errors.

""" Run Phase I """
"-----------------------------------------------------------------------------"
from Phase_1 import phase_1
initial_solution, pass_p2 = phase_1(user_inputs, conversion)

""" Run Phase II """
"-----------------------------------------------------------------------------"
from Phase_2 import phase_2
solution, pass_rhs = phase_2(user_inputs, conversion, initial_solution, pass_p2)

""" Display Requested Variables """
"-----------------------------------------------------------------------------"
if solution['bounded'] == 'yes':
    for i, e in enumerate(var_names):
        print('  ', e, solution['variables'][v.get(e)])

""" If specified, perform a RHS Sensitivity Analysis """
"-----------------------------------------------------------------------------"
if all([sensitivity == 'on', solution['bounded'] == 'yes']):
    from RHS_Sensitivity import rhs_sensitivity
    rhs_sensitivity_report = rhs_sensitivity(user_inputs, constraint, 
                                             conversion, solution, pass_rhs)