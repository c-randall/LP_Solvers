"""
Created on Tue Apr 30 17:18:43 2019

Author: 
    Corey R. Randall

Summary:
    This file pulls in the appropriate user inputs from the LP_Input.py file
    and executes the chosen algorithms with specified options to solve the 
    provided problem.
    
"""

import sys, os
sys.tracebacklimit = 1000 # Set error reporting to default incase different
os.chdir(cwd + '/Simplex_Files')

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
user_inputs['problem'] = problem
user_inputs['dec'] = decimals

from Standard_Form import standard_form
conversion = standard_form(user_inputs, constraint)

if problem == 'dual': # Create and solve the dual problem if specified
    from Dual_Problem import dual_problem
    conversion = dual_problem(user_inputs, conversion)

""" Pre-Processing """
"-----------------------------------------------------------------------------"
# Fill this in with subfunction that simplifies and checks for errors.

""" Import Appropriate Algorithm Subfunctions """
"-----------------------------------------------------------------------------"
if method == 'primal':
    from Phase_1_Primal_Simplex import phase_1
    from Phase_2_Primal_Simplex import phase_2
elif method == 'dual':
    from Phase_1_Dual_Simplex import phase_1
    from Phase_2_Dual_Simplex import phase_2
    
""" Run Phase I and II """
"-----------------------------------------------------------------------------"
initial_solution, pass_p2 = phase_1(user_inputs, conversion)
solution, pass_rhs = phase_2(user_inputs, conversion, initial_solution, pass_p2)

""" Display Requested Variables """
"-----------------------------------------------------------------------------"
if all([solution['bounded'] == 'yes', solution['feasibility'] == 'feasible']):
    if len(var_names) != 0:
        for i, e in enumerate(var_names):
            print('  ', e, np.round(solution['variables'][v.get(e)], decimals))
    else:
        print('   x =', np.round(solution['variables'], decimals))

""" If specified, perform a RHS Sensitivity Analysis """
"-----------------------------------------------------------------------------"
if all([sensitivity == 'on', solution['bounded'] == 'yes', solution['feasibility'] 
        == 'feasible']):
    
    from RHS_Sensitivity import rhs_sensitivity
    rhs_sensitivity_report = rhs_sensitivity(user_inputs, constraint, 
                                             conversion, solution, pass_rhs)