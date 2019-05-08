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
import time, sys
import numpy as np
import scipy.sparse as sp
from math import factorial as fact

""" Function definition """
"-----------------------------------------------------------------------------"
def phase_2(user_inputs, conversion, initial_solution, pass_p2):    
    """ Pre-Processing """
    "-------------------------------------------------------------------------"
    t1 = time.time()    
    
    # Extract dictionary:
    objective = user_inputs['objective']
    tolerance = user_inputs['tolerance']
    problem = user_inputs['problem']
    
    m = conversion['m']
    n = conversion['n']
    n_slack = conversion['n_slack']
    
    A = sp.csr_matrix(conversion['A'])
    b = sp.csr_matrix(conversion['b'])
    c_coeff = sp.csr_matrix(conversion['c_coeff'])
    
    Basis = np.copy(initial_solution['Basis'])
    A_B_inv = sp.csr_matrix.copy(pass_p2['A_B_inv'])
            
    # Determine maximum iterations:
    count = 0
    R, C = A.shape
    max_step = fact(C) // (fact(R)*fact(C - R))  
    
    # If feasible... assume bounded:
    if initial_solution['feasibility'] == 'infeasible':
        bounded = 'n/a'
        start = 0
        count = 0
    else:
        bounded = 'yes'
        start = 1
        count = 0
        c_B = sp.csr_matrix(c_coeff[0, Basis -1])
                   
    # Set the index of the incoming variable and pricing scheme:
    if user_inputs['pricing'] == 'most':
        most = 1
        bland = 0
    elif user_inputs['pricing'] == 'bland':
        most = 0
        bland = 1
        
    if user_inputs['incoming'] == 'first':
        incoming_ind = 0
    elif user_inputs['incoming'] == 'last':
        incoming_ind = -1
        
    """ Simplex Method Loop """
    "-------------------------------------------------------------------------"
    while all([count < max_step, start == 1]):           
        """ Step 1.) Compute current basic feasible solution """
        # Calculate the solution for the current basis:
        b_bar = A_B_inv.dot(b)
        
        # Calculate the dual varialbes and reduced cost:
        y_bar = sp.lil_matrix(c_B.dot(A_B_inv))
        c_bar = sp.lil_matrix(c_coeff - y_bar.dot(A))
        
        """ Step 2.) Check optimality - if needed choose incoming varialbe """                
        if sp.find(c_bar < -1*tolerance)[1].size == 0:
            # If no improvement, print solution to variables and objective:
            ind_Basis = list(Basis -1)
            x = np.zeros([c_coeff.shape[1], c_coeff.shape[0]])
            x[ind_Basis] = b_bar.A
            optimal_value = np.round(c_coeff.dot(x), user_inputs['dec'])
            
            if all([objective == 'maximize', problem == 'primal']):
                optimal_value = -1*optimal_value
            elif all([objective == 'minimize', problem == 'dual']):
                optimal_value = -1*optimal_value
            
            print('STOP (P2): c_bar values were all >= 0. \n')
            
            print('Solution:')
            print('   Optimal Objective Value:', optimal_value[0,0])
            print('   iterations =',count, ', time =',round(time.time()-t1,2), 's')
            
            break
        
        else:
            t_ind = most*np.argmin(c_bar.A)\
                  + bland*sp.find(c_bar < -1*tolerance)[1][incoming_ind]
            
        """ Step 3.) Select outgoing variable via minimum ratio test """
        A_t = A[:, t_ind]
        A_bar_t = A_B_inv.dot(A_t)
                   
        if sp.find(A_bar_t > tolerance)[1].size == 0:
            bounded = 'no'
            print('STOP (P2): problem is unbounded')
            print('iterations =', count, 
                  ', time =', round(time.time()-t1,2), 's')
            break
    
        i_pos = sp.find(A_bar_t > tolerance)[0] # finds A_bar_t > 0            
        ratio_test = b_bar / A_bar_t
        r_ind = i_pos[np.argmin(ratio_test[i_pos])]
        
        """ Step 4.) Update the basis for algorithm to repeat """
        Basis[r_ind] = t_ind +1
        P = sp.identity(m, format='lil')
        P[:,r_ind] = -A_bar_t / A_bar_t[r_ind,0]
        P[r_ind,r_ind] = 1 / A_bar_t[r_ind,0]
        A_B_inv = P.dot(A_B_inv)
        
        c_B = sp.csr_matrix.copy(c_coeff[0, Basis -1])
        
        count = count + 1
        
        if user_inputs['verbose'] == 'on':
            print('P2 count:', count, 'obj:', c_B.dot(b_bar).A[0])
    
    #    user_in = input('"Enter" to continue or "Ctrl+d" to cancel.')   
    #    if user_in == KeyboardInterrupt:
    #        sys.exit(0)
                
    """ Post-Processing for Output Dictionary """
    "-------------------------------------------------------------------------"
    solution = {}
    solution['feasibility'] = initial_solution['feasibility']
    solution['bounded'] = bounded
    solution['count'] = count
    solution['time'] = round(time.time()-t1,2)
    
    pass_rhs = {}
    pass_rhs['A_B_inv'] = A_B_inv
    
    # Store Phase 2 solution if bounded:
    if all([bounded == 'yes', problem == 'primal']):
        if objective == 'maximize':
            y_bar = -1*y_bar 
            
        solution['Objective'] = optimal_value
        solution['Basis'] = Basis
        solution['variables'] = x[:n,0].T
        solution['slacks'] = x[n:,0].T
        solution['duals'] = y_bar.A
        
    elif all([bounded == 'yes', problem == 'dual']): 
        if objective == 'minimize':
            y_bar = -1*y_bar
            
        solution['Objective'] = optimal_value
        solution['Basis'] = Basis
        solution['variables'] = abs(y_bar.A[0,:conversion['n_prim']].T)
        solution['slacks'] = abs(y_bar.A[0,conversion['n_prim']:].T)
        solution['duals'] = x[0:n:2].T - x[1:n:2].T
        
    return solution, pass_rhs
