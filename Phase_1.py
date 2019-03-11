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
import numpy as np
import time

def phase_1(user_inputs, conversion):
    """ Pre-Processing """
    "-------------------------------------------------------------------------" 
    t1 = time.time()
    
    # Extract dictionary for readability:
    A = conversion['A']
    b = conversion['b']
    
    n = conversion['n']
    m = conversion['m']
    n_slack = conversion['n_slack']
    
    # Augment matrix with "artificial variables" as necessary:
    Basis = np.arange(m)
    I = np.identity(m)
    a_ind = n+n_slack
    A_aug = A
    
    for i in range(m):
        
        i_B = []
    
        for j in range(len(A[0,:])):
            if np.sum(abs(I[:,i] - A[:,j])) == 0:
                i_B.append(j +1)
    
        if not i_B:
            i_B.append(a_ind +1)
            a_ind += 1
            A_aug = np.hstack([A_aug, np.reshape(I[:,i],[m,1])])
            
        Basis[i] = i_B[0]

    c_vars = np.zeros(n+n_slack)
    c_arts = np.ones(len(A_aug[0,:]) - len(A[0,:]))
    c_tilde = np.hstack((c_vars, c_arts))
    A_B_inv = I
    
    # Force Loop to Start and Assume Feasible:
    feasibility = 'feasible'
    start = 1
    count = 0
    
    # Build the Basis c-vector:
    c_B = np.zeros(m)
    for i in range(m):
        c_B[i] = c_tilde[Basis[i]-1]
        
    # Set the index of the incoming variable to be first or last neg c_bar:
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
        c_bar = c_tilde - y_bar.dot(A_aug)
        
        """ Step 2.) Check optimality - if needed choose incoming varialbe """
        if all(c_bar >= 0):
            # If no improvement, stop Phase 1 and check feasibility:
            ind_Basis = [int(i)-1 for i in Basis]
            x = np.zeros(len(c_tilde))
            x[ind_Basis] = b_bar[:,0]
            
            if any(c_B != 0):
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
            neg_ind = [i for i, e in enumerate(c_bar) if e < 0]
            t_ind = neg_ind[incoming_ind] + 1
            
        """ Step 3.) Select outgoing variable via minimum ratio test """
        A_t = A_aug[:, t_ind-1]
        A_bar_t = A_B_inv.dot(A_t)
                        
        A_bar_t[A_bar_t == 0] = -1            # Make any 0's negative (ignored)
        ratio_test = b_bar[:,0] / A_bar_t       
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
            c_B[i] = c_tilde[Basis[i]-1]
        
        count = count + 1
    
    """ Post-Processing for Outputs """
    "-------------------------------------------------------------------------"   
    # Generate dictionary for outputs:
    initial_solution = {}
    initial_solution['Basis'] = Basis
    initial_solution['variables'] = x[:n]
    initial_solution['slacks'] = x[n:n+n_slack]
    initial_solution['auxilary'] = x[n+n_slack:]
    initial_solution['feasibility'] = feasibility
    initial_solution['count'] = count       
            
    return initial_solution
