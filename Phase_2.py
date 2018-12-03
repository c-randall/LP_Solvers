"""
Created on Mon Oct 29 16:00:14 2018

Author: 
    Corey R. Randall

Summary:
    Phase II of Simplex Method coded in Python for MEGN 586. Work completed at 
    Colorado School of Mines during the Fall 2018 semester.

"""

""" Import needed modules """
"-----------------------------------------------------------------------------"
import numpy as np

""" Function definition """
"-----------------------------------------------------------------------------"
def phase_2(user_inputs, conversion, initial_solution):
    """ Pre-Processing """
    "-------------------------------------------------------------------------"    
    # Extract dictionary:
    objective = user_inputs['objective']
    tolerance = user_inputs['tolerance']
    
    n = conversion['n']
    m = conversion['m']
    
    A = conversion['A']
    b = conversion['b']
    c_coeff = conversion['c_coeff']
    
    Basis = initial_solution['Basis']
    A_B_inv = initial_solution['A_B_inv']  
    
    # If Feasible... Force Loop to Start:
    if initial_solution['feasibility'] == 'infeasible':
        start = 0
    else:
        start = 1
        count = 1
    
    # Build the Basis c-vector:
        c_B = np.zeros(m)
        for i in range(m):
            c_B[i] = c_coeff[Basis[i]-1]
            
            
    # Set the index of the incoming variable to be first or last neg c_bar:
    if user_inputs['incoming'] == 'first':
        incoming_ind = 0
    elif user_inputs['incoming'] == 'last':
        incoming_ind = -1
    
    P = np.identity(m)
    
    """ Simplex Method Loop """
    "-------------------------------------------------------------------------"
    while start == 1:       
        """ Step 1.) Compute current basic feasible solution """
        # Calculate the solution for the current basis:
        b_bar = A_B_inv.dot(b)
        
        # Calculate the dual varialbes and reduced cost:
        y_bar = c_B.dot(A_B_inv)
        c_bar = c_coeff - y_bar.dot(A)
        
        """ Step 2.) Check optimality - if needed choose incoming varialbe """
        c_bar[abs(c_bar) <= tolerance] = 0
        
        if all(c_bar >= 0):
            # If no improvement, print solution to variables and objective:
            ind_Basis = [int(i)-1 for i in Basis]
            x = np.zeros(len(c_coeff))
            x[ind_Basis] = b_bar[:,0]
            optimal_value = np.round(sum(c_coeff*x), user_inputs['dec'])
            
            if objective == 'maximize':
                optimal_value = -1*optimal_value
            
            print('STOP (P2): c_bar values were all >= 0. \n')
            print('Solution: \n x =', np.round(x[0:n], user_inputs['dec']))
            print('\n Optimal Objective Value:', optimal_value)
            break
        
        else:
            neg_ind = [i for i, e in enumerate(c_bar) if e < 0]
            t_ind = neg_ind[incoming_ind] + 1
            
        """ Step 3.) Select outgoing variable via minimum ratio test """
        A_t = A[:, t_ind-1]
        A_bar_t = A_B_inv.dot(A_t)
        A_bar_t_ratio = A_bar_t
                   
        if all(A_bar_t <= 0):
            print('STOP (P2): problem is unbounded')
            break
    
        A_bar_t_ratio[A_bar_t == 0] = -1      # Make any 0's negative (ignored)
        ratio_test = b_bar[:,0] / A_bar_t_ratio
        min_ratio = min([n for n in ratio_test if n >= 0])
        r_ind = [i for i, e in enumerate(ratio_test) if e == min_ratio][0]
        
        A_bar_t = A_B_inv.dot(A_t)  # Ensure no overwrites to A_bar_t were made
        
        """ Step 4.) Update the basis for algorithm to repeat """
        Basis[r_ind] = t_ind
        P = np.identity(m)
        P[:,r_ind] = -A_bar_t / A_bar_t[r_ind]
        P[r_ind,r_ind] = 1 / A_bar_t[r_ind]
        A_B_inv = P.dot(A_B_inv)
        
        for i in range(m):
            c_B[i] = c_coeff[Basis[i]-1]
        
        count = count + 1
                
    """ Post-Processing for Output Dictionary """
    "-------------------------------------------------------------------------"
    solution = {}
    if objective == 'maximize':
        y_bar = -1*y_bar
    
    # If Feasible... Store Phase 2 Solution:
    if initial_solution['feasibility'] == 'feasible':
        solution['A_B_inv'] = A_B_inv
        solution['duals'] = y_bar
        solution['variables'] = x[:n]
        solution['slacks'] = x[n:]
        solution['Basis'] = Basis
        solution['count_P2'] = count
    
    return solution