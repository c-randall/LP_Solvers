"""
Created on Tue Oct 30 09:18:34 2018

Author: 
    Corey R. Randall

Summary:
    Standard Form conversion tool written for MEGN 586. This function takes the
    information from the LP_Runner file and converts the user input expressions
    into standard form for the Simplex Method to solve. Work completed at 
    Colorado School of Mines during the Fall 2018 semester.

"""

""" Import needed modules """
"-----------------------------------------------------------------------------"
import numpy as np

""" Function definition """
"-----------------------------------------------------------------------------"
def standard_form(user_inputs, constraint):
    """ Pre-Processing """
    "-------------------------------------------------------------------------"    
    # Extract dictionary: 
    objective = user_inputs['objective']
    c_coeff = np.array(user_inputs['c_coeff'])
    pointer = user_inputs['pointer']
    
    # Determine n and m:
    n = len(c_coeff)        # number of variables (excluding slack/excess)
    m = len(constraint)     # number of constraints (excluding non-negativity)
    
    n_slack = 0             # number of slack/excess variables
    for i in range(m):
        if constraint[i+1][pointer['inequality']] == '<=':
            n_slack = n_slack + 1
        elif constraint[i+1][pointer['inequality']] == '>=':
            n_slack = n_slack + 1
    
    # Initialize size of A matrix:
    A = np.zeros([m, n+n_slack])
    
    """ Construct and Convert Problem """
    "-------------------------------------------------------------------------"
    # Convert to 'minimize' for standard form:
    if objective == 'maximize':
        c_coeff = -1*c_coeff
    
    if n_slack != 0:
        c_coeff = np.hstack((c_coeff, np.zeros(n_slack)))
    
    # Construct A Matrix and b vector:
    b = np.zeros([m, 1])
    for i in range(m):
        extract = constraint[i+1]
        
        if extract[pointer['inequality']] == '<=':
            A[i,n+i] = 1
        elif extract[pointer['inequality']] == '>=':
            A[i,n+i] = -1
    
    # Ensure that equations have positive b-values:
        A[i,0:n] = np.array(extract[pointer['coeff']])
        b[i,0] = extract[pointer['b-value']]
        if extract[pointer['b-value']] < 0:
            A[i,:] = -1*A[i,:]
            b[i,0] = -1*b[i,0]         
                
    """ Post-Processing for Outputs """
    "-------------------------------------------------------------------------"    
    # Generate dictionary for outputs:
    conversion = {}
    conversion['A'] = A
    conversion['b'] = b
    conversion['c_coeff'] = np.reshape(c_coeff, [1, c_coeff.size])
    
    conversion['n'] = n
    conversion['m'] = m
    conversion['n_slack'] = n_slack
    
    return conversion
