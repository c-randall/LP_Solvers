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

def phase_1(user_inputs, conversion):
    """ Pre-Processing """
    "-------------------------------------------------------------------------"   
    # Extract dictionary for readability:
    A = conversion['A']
    b = conversion['b']
    
    n = conversion['n']
    m = conversion['m']
    n_slack = conversion['n_slack']
    
    # Augment matrix with "artificial variables" as starting Basis:
    A_aug = np.hstack((A, np.identity(m)))
    c_tilde = np.hstack((np.zeros(n+n_slack), np.ones(m)))
    Basis = np.arange(n+n_slack, n+n_slack+m) + 1
    A_B_inv = np.identity(m)
    
    # Force Loop to Start and Assume Feasible:
    feasibility = 'feasible'
    start = 1
    count = 1
    
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
    # Erase Later:
    ###########################################################################
#    Basis = [1, 2, 5]
#    A_B_inv = np.array([[1/2, -1/2, 0], [1/2, 1/2, 0], [-1/2, -1/2, 1]])
#    feasibility = 'infeasible'
    ###########################################################################
    
    # Generate dictionary for outputs:
    initial_solution = {}
    initial_solution['Basis'] = Basis
    initial_solution['A_B_inv'] = A_B_inv
    initial_solution['solution'] = x
    initial_solution['variables'] = x[:n]
    initial_solution['slacks'] = x[n:n+n_slack]
    initial_solution['auxilary'] = x[n+n_slack:]
    initial_solution['feasibility'] = feasibility
    initial_solution['count_P1'] = count
    
    if feasibility == 'feasible':
        print('End Phase 1... begin Phase 2. \n')
    
    return initial_solution
