"""
Created on Wed Mar 20 13:34:14 2019

Author: 
    Corey R. Randall

Summary:
    Phase II of Dual Simplex Method coded in Python.

"""

""" Import needed modules """
"-----------------------------------------------------------------------------"
import scipy.sparse as sp
import numpy as np
import time, warnings

warnings.simplefilter('ignore', sp.SparseEfficiencyWarning)

""" Function definition """
"-----------------------------------------------------------------------------"
def phase_2(user_inputs, conversion, initial_solution, pass_p2):    
    """ Pre-Processing """
    "-------------------------------------------------------------------------"
    t1 = time.time()    
    
    # Extract dictionary:
    objective = user_inputs['objective']
    tolerance = user_inputs['tolerance']
    
    n = conversion['n']
    m = conversion['m']
    n_slack = conversion['n_slack']
    
    A = sp.csr_matrix(conversion['A'])
    b = sp.csr_matrix(conversion['b'])
    c_coeff = sp.csr_matrix(conversion['c_coeff'])
    
    Basis = np.copy(initial_solution['Basis'])
    A_B_inv = pass_p2['A_B_inv']
    
    # If Feasible... Force Loop to Start:
    if initial_solution['bounded'] == 'no':
        start = 0
        count = 0
        feasibility = 'infeasible'
    else:
        start = 1
        count = 0
        feasibility = 'feasible'
    
    # Build the Basis c-vector:
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
    while start == 1:
        """ Step 1.) Compute current basic feasible solution """
        # Calculate the solution for the current basis:
        b_bar = A_B_inv.dot(b)
        b_bar[np.where(abs(b_bar.A) < tolerance)[0], 0] = 0
        
        # Calculate the dual varialbes and reduced cost:
        y_bar = c_B.dot(A_B_inv)
        c_bar = c_coeff - y_bar.dot(A)
        c_bar[0, np.where(abs(c_bar.A) < tolerance)[1]] = 0
        
        """ Step 2.) Check optimality - if needed choose outgoing varialbe """        
        if sp.find(b_bar < 0)[0].size == 0:
            # If no improvement, print solution to variables and objective:
            ind_Basis = list(Basis -1)
            x = np.zeros([c_coeff.shape[1], c_coeff.shape[0]])
            x[ind_Basis] = b_bar.A
            optimal_value = np.round(c_coeff.dot(x), user_inputs['dec'])
            
            if objective == 'maximize':
                optimal_value = -1*optimal_value
            
            print('STOP (P2): b_bar values were all >= 0. \n')
            
            print('Solution:')
            print('   Optimal Objective Value:', optimal_value[0,0])
            print('   iterations =',count, ', time =',round(time.time()-t1,2), 's')
            
            break
        
        else:
            r_ind = most*np.argmin(b_bar.A) + bland*sp.find(b_bar < 0)[0][incoming_ind]
            
        """ Step 3.) Select incoming variable via minimum ratio test """
        A_bar_r = A_B_inv[r_ind,:].dot(A)
        
        if sp.find(A_bar_r < 0)[1].size == 0:
            feasibility = 'infeasible'
            print('STOP (P2):', feasibility, '- dual is unbounded')
            print('iterations =', count, 
                  ', time =', round(time.time()-t1,2), 's')
            break
    
        i_neg = sp.find(A_bar_r < -tolerance)[1] # finds A_bar_t < 0
        ratio_test = -c_bar / A_bar_r
        t_ind = i_neg[np.argmin(ratio_test[0,i_neg])]                 
        
        """ Step 4.) Update the basis for algorithm to repeat """
        A_t = A[:, t_ind]
        A_bar_t = A_B_inv.dot(A_t)
        
        Basis[r_ind] = t_ind +1
        P = sp.identity(m, format='lil')
        P[:,r_ind] = -A_bar_t / A_bar_t[r_ind,0]
        P[r_ind,r_ind] = 1 / A_bar_t[r_ind,0]
        A_B_inv = P.dot(A_B_inv)
        
        c_B = c_coeff[0, Basis -1]
        
        count = count + 1
        
    #    print('P2 count:', count, 'obj:', c_B.dot(b_bar).A[0])
    
    #    user_in = input('"Enter" to continue or "Ctrl+d" to cancel.')   
    #    if user_in == KeyboardInterrupt:
    #        sys.exit(0)
                
    """ Post-Processing for Output Dictionary """
    "-------------------------------------------------------------------------"
    solution = {}
    solution['feasibility'] = feasibility
    solution['bounded'] = initial_solution['bounded']
    solution['count'] = count
    solution['time'] = round(time.time()-t1,2)
    
    pass_rhs = {}
    pass_rhs['A_B_inv'] = A_B_inv
    
    # Store Phase 2 solution if feasible:
    if feasibility == 'feasible':
        
        solution['Objective'] = optimal_value
        solution['Basis'] = Basis
        if objective == 'maximize':
            y_bar = -1*y_bar
        
        solution['variables'] = x[:n,0].T
        solution['slacks'] = x[n:n+n_slack,0].T
        solution['excess'] = x[n+n_slack:,0].T
        solution['duals'] = y_bar.A
        
    return solution, pass_rhs