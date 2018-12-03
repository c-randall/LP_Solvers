"""
Created on Wed Nov 28 12:12:23 2018

Author: 
    Corey R. Randall

Summary:
    RHS sensitivity analysis for LPs coded in Python for MEGN 586. Work 
    completed at Colorado School of Mines during the Fall 2018 semester.

"""

""" Import needed modules """
"-----------------------------------------------------------------------------"
import numpy as np

""" Function definition """
"-----------------------------------------------------------------------------"
def rhs_sensitivity(user_inputs, constraint, conversion, solution):
    
    # Read out user values and extract current optimal solution:
    pointer = user_inputs['pointer']
    Nc = len(constraint)
    
    x_all = np.hstack([solution['variables'], solution['slacks']])
    x_var = solution['variables']
    b = conversion['b']
    
    # Loop over each constraint to change its respective rhs:
    Non_bind_ind = [ind for ind, e in enumerate(solution['slacks']) if e != 0]
    
    rhs_sensitivity_report = {}
    
    for i in range(Nc):
        e_i = np.zeros([Nc,1])
        e_i[i,0] = 1.
        
        # Figure out how basis varialbes change:
        x_change = solution['A_B_inv'].dot(e_i)
        x_order = np.concatenate((solution['Basis'].reshape(Nc,1), 
                                  x_change), axis=1)
        x_order.sort(axis=0)

        # Pull variable changes (neglecting slacks) for constraint checks:
        x_change_ind = [ind for ind, e in enumerate(x_order[:,0]) if e <= Nc]
        x_change_vars = np.zeros_like(solution['variables'])
        x_change_vars[x_change_ind] = x_order[x_change_ind,1]
        
        # Initialize range vectors (will be sorted later):
        delt_sign = []
        delt_val = np.zeros(len(Non_bind_ind)+len(solution['Basis']))
        
        # Loop over non-binding constraints to maintain feasibility:
        count = 0
        
        for j, NB_ind in enumerate(Non_bind_ind):
            coeff_const = constraint[NB_ind+1][pointer['coeff']]
            
            delta_coeff = np.sum(coeff_const*x_change_vars)
            
            if delta_coeff != 0:
                if delta_coeff < 0.:
                    if constraint[NB_ind+1][pointer['inequality']] == '>=':
                        delt_sign.append('<=')
                    elif constraint[NB_ind+1][pointer['inequality']] == '<=':
                        delt_sign.append('>=')
                else:
                    delt_sign.append(
                              constraint[NB_ind+1][pointer['inequality']])
                    
                delt_val[count] = (b[NB_ind] - np.sum(coeff_const*x_var))\
                                / delta_coeff
                                
                count = count +1
            
        # Loop over changed variables to maintain non-negativity:
        for j, ind_B in enumerate(solution['Basis']):
            if x_change[j] != 0:
                if x_change[j] < 0:
                    delt_sign.append('<=')
                else:
                    delt_sign.append('>=')
                    
                delt_val[count] = -x_all[ind_B-1] / x_change[j]
                count = count +1
                
        # Determine most restrictive upper and lower bounds:
        LB_ind = [ind for ind, e in enumerate(delt_sign) if e == '>=']
        UB_ind = [ind for ind, e in enumerate(delt_sign) if e == '<=']
        
        # Report infinite values if no upper or lower bounds:
        if len(LB_ind) == 0:
            LB = '-inf'
            rhs_down = '-inf'
        else:
            LB = np.max(delt_val[LB_ind])
            rhs_down = np.round(float(b[i] + LB), user_inputs['dec'])
            
        if len(UB_ind) == 0:
            UB = 'inf'
            rhs_up = 'inf'
        else:
            UB = np.min(delt_val[UB_ind])
            rhs_up = np.round(float(b[i] + UB), user_inputs['dec'])
        
        rhs_sensitivity_report[i] = constraint[i+1][pointer['name']]+'_range',\
                            [rhs_down, rhs_up], 'dual:', \
                            np.round(solution['duals'][i], user_inputs['dec'])
    
    return rhs_sensitivity_report