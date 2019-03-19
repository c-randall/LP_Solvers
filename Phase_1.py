"""
Created on Mon Oct 29 15:27:08 2018

Author: 
    Corey R. Randall

Summary:
    Phase I of Simplex Method coded in Python for MEGN 586. Work completed at 
    Colorado School of Mines during the Fall 2018 semester.

"""

""" Import needed modules """
"-----------------------------------------------------------------------------"
import scipy.sparse as sp
import numpy as np
import time, sys

def phase_1(user_inputs, conversion):
    """ Pre-Processing """
    "-------------------------------------------------------------------------" 
    t1 = time.time()
    
    # Extract dictionary for readability:
    tolerance = user_inputs['tolerance']
    
    n = conversion['n']
    m = conversion['m']
    n_slack = conversion['n_slack']
    
    A = sp.csc_matrix(conversion['A'])
    b = sp.csc_matrix(conversion['b'])
    
    # Augment matrix with "artificial variables" as necessary:  
    I_ind = []
    i_B = []
    
    " Locate columns that have a single positive 1 - store column and row "
    for i in range(n + n_slack):
        col = A[:,i]
        
        if all([abs(sp.csc_matrix.sum(col)) == 1, sp.csc_matrix.sum(col) == 1]):
            if sp.find(col == 1)[0] not in I_ind:
                I_ind.append(int(sp.find(col == 1)[0]))
                i_B.append(i)
    
    " Insert artificials for any rows missing 1 from an identity column "     
    i_art = list(set(range(m)) - set(I_ind))
    i_B.extend(range(n + n_slack, n + n_slack + len(i_art))) 
    
    " Build columns to tack onto A matrix from arificial variables "
    A_art = np.zeros([m, len(i_art)])
    for i, e in enumerate(i_art):
        A_art[e,i] = 1

    " Start with identity inverse and augment A and c with artificials "       
    A_B_inv = sp.identity(m)    
    A_aug = sp.hstack([A, A_art], format='csc')
    c_tilde = sp.hstack([np.zeros(n+n_slack), np.ones(len(i_art))], format='csc')
    
    " Arrange basis to match locations of identity columns "
    I_ind.extend(i_art)
    Basis = np.vstack([I_ind, np.array(i_B) +1])
    Basis = Basis[-1, Basis[0].argsort()]
    
    # Force Loop to Start and Assume Feasible:
    feasibility = 'feasible'
    start = 1
    count = 0
    
    # Build the Basis c-vector:
    c_B = sp.csr_matrix.copy(c_tilde[0, Basis -1])
        
    # Set the index of the incoming variable for bland's rule:        
    if user_inputs['incoming'] == 'first':
        incoming_ind = 0
    elif user_inputs['incoming'] == 'last':
        incoming_ind = -1
            
    """ Begin Simplex Loop to Replace Artificial Variables """
    "-------------------------------------------------------------------------"
    while start == 1:         
        """ Step 1.) Compute current basic feasible solution """
        # Calculate the solution for the current basis:
        b_bar = A_B_inv.dot(b)
        
        # Calculate the dual varialbes and reduced cost:
        y_bar = c_B.dot(A_B_inv)
        c_bar = sp.lil_matrix(c_tilde - y_bar.dot(A_aug))
        
        """ Step 2.) Check optimality - if needed choose incoming varialbe """
        c_bar[0, np.where(abs(c_bar.A) < tolerance)[1]] = 0
    
        if sp.find(c_bar < 0)[1].size == 0:
            # If no improvement, stop Phase 1 and check feasibility:
            ind_Basis = list(Basis -1)
            x = np.zeros([c_tilde.shape[1], c_tilde.shape[0]])
            x[ind_Basis] = b_bar.A
            
            if sp.find(c_B > 0)[1].size != 0:
                feasibility = 'infeasible'
                print('STOP (P1):', feasibility, '- artificials in Basis')
                print('iterations =', count, 
                      ', time =', round(time.time()-t1,2), 's \n')
                
            else:
                print('End Phase 1... begin Phase 2.')
                print('iterations =', count, 
                      ', time =', round(time.time()-t1,2), 's \n')
            
            break
        
        else:
            t_ind = sp.find(c_bar < 0)[1][incoming_ind] + 1
            
        """ Step 3.) Select outgoing variable via minimum ratio test """
        A_t = A_aug[:, t_ind -1]
        A_bar_t = A_B_inv.dot(A_t)
                        
        i_pos = sp.find(A_bar_t > tolerance)[0] # finds A_bar_t > 0
        ratio_test = -1*np.ones(b_bar.shape)            
        ratio_test[i_pos] = b_bar[i_pos,0] / A_bar_t[i_pos,0]
        min_ratio = min([n for n in ratio_test if n >= 0])
        r_ind = np.where(ratio_test == min_ratio)[0][0]
        
        """ Step 4.) Update the basis for algorithm to repeat """
        Basis[r_ind] = t_ind
        P = sp.identity(m, format='lil')
        P[:,r_ind] = -A_bar_t / A_bar_t[r_ind,0]
        P[r_ind,r_ind] = 1 / A_bar_t[r_ind,0]
        A_B_inv = P.dot(A_B_inv)
        
        c_B = sp.csr_matrix.copy(c_tilde[0, Basis -1])
        
        count = count + 1
    
    #    print('P1 count:', count, 'obj:', c_B.dot(b_bar).A[0]) 
    
    #    user_in = input('"Enter" to continue or "Ctrl+d" to cancel.')   
    #    if user_in == KeyboardInterrupt:
    #        sys.exit(0)       
            
    """ Post-Processing for Outputs """
    "-------------------------------------------------------------------------"   
    # Generate dictionary for outputs:
    initial_solution = {}
    initial_solution['Basis'] = Basis
    initial_solution['variables'] = x[:n,0].T
    initial_solution['slacks'] = x[n:n+n_slack,0].T
    initial_solution['feasibility'] = feasibility
    initial_solution['count'] = count
    initial_solution['time'] = round(time.time()-t1,2)      
    
    pass_p2 = {}
    pass_p2['A_B_inv'] = A_B_inv
            
    return initial_solution, pass_p2
