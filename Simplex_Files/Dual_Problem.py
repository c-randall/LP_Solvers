"""
Created on Wed Apr  3 13:07:18 2019

Author: 
    Corey R. Randall

Summary:
    If the user wishes to solve the Dual Problem over the Primal one, this 
    function provides support to appropriately convert the problem into its
    alternate form.

"""

""" Import needed modules """
"-----------------------------------------------------------------------------"
import numpy as np

""" Function definition """
"-----------------------------------------------------------------------------"
def dual_problem(user_inputs, conversion):
    # Extract dictionary for readibility:
    A = conversion['A']
    b = conversion['b']
    c_coeff = conversion['c_coeff']
    
    n = conversion['n']
    m = conversion['m']
    n_slack = conversion['n_slack']
    
    # Convert A, c_coeff to allow for unrestricted y (i.e. y = y' - y''):
    A_temp = np.repeat(A.T, 2)
    A_temp[1::2] = -A_temp[1::2]
    A = np.reshape(A_temp, [A.shape[1], 2*A.shape[0]])
    A = np.hstack([A, np.identity(A.shape[0])])
    
    b_temp = np.repeat(b.T, 2) # the obj. coeff. in (D) are b from (P)
    b_temp[1::2] = -b_temp[1::2]
    b = np.reshape(b_temp, [b.shape[1], 2*b.shape[0]])
    b = np.hstack([b, np.zeros([1, A.shape[0]])])
    
    # Ensure no negative values on RHS:
    for i in range(c_coeff.shape[1]):
        if c_coeff[0,i] < 0: # the RHS, b values, in (D) are c from (P) 
            A[i,:] = -A[i,:]
            c_coeff[0,i] = -c_coeff[0,i]
    
    # Generate dictionary for outputs:
    dual_conversion = {}
    dual_conversion['A'] = A
    dual_conversion['b'] = c_coeff.T
    dual_conversion['c_coeff'] = -b
    
    dual_conversion['n'] = 2*m
    dual_conversion['m'] = n +n_slack
    dual_conversion['n_slack'] = n +n_slack
    
    dual_conversion['n_prim'] = n
    dual_conversion['n_slack_prim'] = n_slack
    
    return dual_conversion
