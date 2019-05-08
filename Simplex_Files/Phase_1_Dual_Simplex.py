"""
Created on Wed Mar 20 13:33:21 2019

Author: 
    Corey R. Randall

Summary:
    Phase I of Simplex Method coded in Python for MEGN 586. Work completed at 
    Colorado School of Mines during the Fall 2018 semester.

"""

""" Import needed modules """
"-----------------------------------------------------------------------------"
import scipy, time, sys
import numpy as np
import scipy.sparse as sp
import scipy.linalg as linalg
from math import factorial as fact

""" Function definition """
"-----------------------------------------------------------------------------"
def phase_1(user_inputs, conversion):
    """ Pre-Processing """
    "-------------------------------------------------------------------------" 
    t1 = time.time()
    
    # Extract dictionary for readability:    
    n = conversion['n']
    m = conversion['m']
    n_slack = conversion['n_slack']
    
    A = sp.csr_matrix(conversion['A'])
    b = sp.csr_matrix(conversion['b'])
    c_coeff = sp.csr_matrix(conversion['c_coeff'])
    
    # Exit if any c_coeff < 0:
    if sp.find(c_coeff < 0)[1].size != 0:
        sys.tracebacklimit = 0
        raise Exception('Use Primal Simplex due to c_coeff structure.')
        sys.tracebacklimit = 1000
        
    # Attempt to find c_B with all 0 coefficients:
    if c_coeff.shape[1] - sp.find(c_coeff != 0)[1].size >= m:
        i_N = sp.find(c_coeff != 0)[1]
        i_B = list(set(range(n +n_slack)) - set(i_N))[:m]
        
        # Build the Basis from the excess variables:
        Basis = np.array(i_B) +1
    
        # Set bounded:
        bounded = 'yes'
        count = 0
    
        # Calculate initial solution values:
        A_B_inv = sp.csr_matrix(np.linalg.inv(A.A[:,i_B]))
        b_bar = A_B_inv.dot(b)
        x = np.zeros([n +n_slack, 1])
        x[Basis -1] = b_bar.A
        
        print('End Phase 1... begin Phase 2.')
        print('iterations =', count, ', time =', round(time.time()-t1,2), 's \n')
    
    # Convert to Dual problem and use Primal Phase I to get starting basis:    
    else:
        from Dual_Problem import dual_problem
        conversion_dual = dual_problem(user_inputs, conversion)
        
        # Extract dictionary for readability:
        objective = user_inputs['objective']
        tolerance = user_inputs['tolerance']
        
        n = conversion_dual['n']
        m = conversion_dual['m']
        n_slack = conversion_dual['n_slack']
        
        A = sp.csr_matrix(conversion_dual['A'])
        b = sp.csr_matrix(conversion_dual['b'])
        
        # Augment matrix with "artificial variables" as necessary:  
        I_ind = []
        i_B = []
        
        " Locate columns that have a single positive 1 - store column and row "
        reg_sum = sp.find(np.sum(A, axis=0) == 1)
        abs_sum = sp.find(np.sum(abs(A), axis=0) == 1)
        col_id = np.intersect1d(reg_sum[1], abs_sum[1])
        row_id = sp.find(A[:,col_id] == 1)
        
        I_ind.extend(row_id[0])
        i_B.extend(col_id)
        
        " Insert artificials for any rows missing 1 from an identity column "     
        i_art = list(set(range(m)) - set(I_ind))
        i_B.extend(range(n + n_slack, n + n_slack + len(i_art))) 
        
        " Build columns to tack onto A matrix from arificial variables "
        A_art = np.zeros([m, len(i_art)])
        try:
            A_art[np.array(i_art), np.arange(0,A_art.shape[1])] = 1
        except:
            A_art = A_art
        
        " Start with identity inverse and augment A and c with artificials "       
        A_B_inv = sp.identity(m)    
        A_aug = sp.hstack([A, A_art], format='csr')
        c_tilde = sp.hstack([np.zeros(n+n_slack), np.ones(len(i_art))], format='csr')
        
        " Arrange basis to match locations of identity columns "
        I_ind.extend(i_art)
        Basis = np.vstack([I_ind, np.array(i_B) +1])
        Basis = Basis[-1, Basis[0].argsort()]
        
        # Force Loop to Start and Assume Feasible:
        bounded = 'yes'
        count = 0
        
        # Build the Basis c-vector:
        c_B = sp.csr_matrix.copy(c_tilde[0, Basis -1])
        
        # Determine maximum iterations:
        R, C = A_aug.shape
        max_step = fact(C) // (fact(R)*fact(C - R))
            
        # Set the index of the incoming variable for bland's rule:        
        if user_inputs['incoming'] == 'first':
            incoming_ind = 0
        elif user_inputs['incoming'] == 'last':
            incoming_ind = -1
              
        """ Begin Simplex Loop to Replace Artificial Variables """
        "---------------------------------------------------------------------"
        for i in range(max_step):  
            """ Step 1.) Compute current basic feasible solution """
            # Calculate the solution for the current basis:
            b_bar = A_B_inv.dot(b)
            
            # Calculate the dual varialbes and reduced cost:
            y_bar = sp.lil_matrix(c_B.dot(A_B_inv))
            c_bar = sp.lil_matrix(c_tilde - y_bar.dot(A_aug))
            
            """ Step 2.) Check optimality - if needed choose incoming varialbe """       
            if sp.find(c_bar < -1*tolerance)[1].size == 0:
                # If no improvement, stop Phase 1 and check feasibility:
                x = np.zeros([c_tilde.shape[1], c_tilde.shape[0]])
                x[Basis -1] = b_bar.A
                
                if sp.find(c_B > 0)[1].size != 0:
                    bounded = 'maybe'
                    feasibility = 'infeasible'
                    print('STOP (P1): Dual is', feasibility)
                    print('iterations =', count, 
                          ', time =', round(time.time()-t1,2), 's \n')
                    
                else:
                    print('End Phase 1... begin Phase 2.')
                    print('iterations =', count, 
                          ', time =', round(time.time()-t1,2), 's \n')
                
                break
            
            else:
                t_ind = sp.find(c_bar < -1*tolerance)[1][incoming_ind]
                
            """ Step 3.) Select outgoing variable via minimum ratio test """
            A_t = A_aug[:, t_ind]
            A_bar_t = A_B_inv.dot(A_t)
                            
            i_pos = sp.find(A_bar_t > tolerance)[0] # finds A_bar_t > 0            
            ratio_test = b_bar / A_bar_t
            min_ind = i_pos[np.where(ratio_test[i_pos] == min(ratio_test[i_pos]))[0]]
            
            try: # attempt to replace any aritificial variables first
                inter = np.intersect1d(Basis[Basis > n+n_slack], Basis[min_ind])[0]
                r_ind = np.where(Basis == inter)[0][0]       
            except:
                r_ind = min_ind[-1]
            
            """ Step 4.) Update the basis for algorithm to repeat """
            Basis[r_ind] = t_ind +1
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
        
        # Convert back to using primal problem structure:
        A = sp.csr_matrix(conversion['A'])
        b = sp.csr_matrix(conversion['b'])
        c_coeff = sp.csr_matrix(conversion['c_coeff'])
        
        n = conversion['n']
        n_slack = conversion['n_slack']
        
        if objective == 'minimize':
            y_bar = -1*y_bar
             
        x = y_bar.A.T
        
        try:
            Basis = np.where(x != 0)[0] +1
        except:
            Basis = []
            
        if len(Basis) < conversion['m']:
            z_col = np.where(x == 0)[0]
            ind_max = np.argmax(abs(A[:,z_col]), axis=0)
            A_temp = A[:,z_col] / A.A[ind_max,z_col][0][np.newaxis,:]
            uniq_col = np.unique(A_temp, return_index=True, axis = 1)[1]
            
            z_Basis = z_col[uniq_col[:conversion['m'] - len(Basis)]]
            Basis = np.hstack([Basis, z_Basis])
        
        # Calculate initial solution values:
        try:            
            A_B_inv = sp.csr_matrix(np.linalg.inv(A.A[:, Basis -1]))
            b_bar = A_B_inv.dot(b)
            x = np.zeros([n +n_slack, 1])
            x[Basis -1] = b_bar.A
        except:
            P, L, U = linalg.lu(A.A[:, Basis -1])
            
            try:
                L_inv = np.linalg.inv(L)
            except:
                L_inv = np.linalg.pinv(L)
                
            try:
                U_inv = np.linalg.inv(U)
            except:
                U_inv = np.linalg.pinv(U)
                
            if all([np.allclose(L, np.dot(L, np.dot(L_inv, L)), atol=user_inputs['tolerance']),
                    np.allclose(L_inv, np.dot(L_inv, np.dot(L, L_inv)), atol=user_inputs['tolerance']),
                    np.allclose(U, np.dot(U, np.dot(U_inv, U)), atol=user_inputs['tolerance']),
                    np.allclose(U_inv, np.dot(U_inv, np.dot(U, U_inv)), atol=user_inputs['tolerance'])]):
    
                A_B_inv = sp.csr_matrix(U_inv.dot(L_inv))
                b_bar = A_B_inv.dot(b)
                x = np.zeros([n +n_slack, 1])
                x[Basis -1] = b_bar.A
                
            else:
                sys.tracebacklimit = 0
                raise Exception('Use Primal Simplex due to singular inversion.')
                
        c_B = sp.csr_matrix(c_coeff[0, Basis -1])
        
        y_bar = c_B.dot(A_B_inv)
        c_bar = c_coeff - y_bar.dot(A)
        c_bar[0, np.where(abs(c_bar.A) < tolerance)[1]] = 0
        
#        if (c_bar.A < -1*tolerance).any():
#            sys.tracebacklimit = 0
#            raise Exception('Use Primal Simplex due initialization.')
        
    """ Post-Processing for Outputs """
    "-------------------------------------------------------------------------"   
    # Generate dictionary for outputs:
    initial_solution = {}
    initial_solution['Basis'] = Basis
    initial_solution['variables'] = x[:n,0].T
    initial_solution['slacks'] = x[n:n+n_slack,0].T
    initial_solution['bounded'] = bounded
    initial_solution['count'] = count
    initial_solution['time'] = round(time.time()-t1,2)      
    
    pass_p2 = {}
    pass_p2['A_B_inv'] = A_B_inv
            
    return initial_solution, pass_p2
