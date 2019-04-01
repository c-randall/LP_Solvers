"""
Created on Fri Mar 28 8:33:02 2019

Author: 
    Corey R. Randall

Summary:
    For larger and more complex problems, users may choose to write a separate 
    script that contains the objective function, variables, and constraints. 
    The formulation given below was developed by Mohammad Fathollahzadeh and 
    Karl Heine at Colorado School of Mines. It attempts to size a PV and 
    battery system for Hill Hall (one of the buildings on the CSM campus). This 
    model and its results were presented at the 2019 ASHRAE Winter Conference 
    that took place in Atlanta, GA. With their permission, the problem was 
    translated to python by Corey R. Randall for performance tests on the 
    Primal Simplex Method Algorithm. The more complex problem prompted the 
    addition of sparse matrices and steep pricing to improve the algorithm's 
    performance. Findings showed that these additions reduced solution time 
    from 10 min to 40 sec. Solutions presented at ASHRAE have also been matched
    by the output of this python model and algorithm for validation purposes.
    
    This is a copy of the pv_batt_sizing.py script with substitutions made for
    the penalty varialbes to reduce the number of equations and constraints.
    Without substitutions a combination of Phase 1 and 2 ran in approximately
    26s. These substitutions should elimate at least 1460 variables and 
    constraints to reduce the problem size and overall computation time. After
    substitutions were made, this in fact provided improved efficiency with
    an overall run time of about 8s. This justifies the addition of a presolve
    at a later point in time that cleans up inefficiencies from user input and
    performs variable substitutions when possible (i.e. RHS = 0).
    
"""

""" Import needed modules """
"-----------------------------------------------------------------------------"
import numpy as np

""" Input model """
"-----------------------------------------------------------------------------"
# Sets:
T = 365

# Parameters:
" Constant parameters - i.e. do not change based on index "
C_B = 385       # batt cost [$/kWh storage] - SAM (5MWh batt, 700kW max power)
C_PV = 1820     # PV cost [$/kW generation] - SAM (14MW PV system)
C_e = 0.066917  # electricity cost [$/kWh] - current Mines cost to Xcel
Y_B = 8         # expected battery life [yr] - previous LiNMC studies
Y_PV = 25       # expected PV life [yr] - based on SAM default
D_s = 16.578    # summer demand charge [$/kW] - current Mines cost to Xcel
D_w = 14.777    # winter demand charge [$/kW] - current Mines cost to Xcel
E_B = 0.85      # battery efficiency [-] - round trip estimate for LiNMC
E_PV = 1.0      # PV efficiency [-] - accounted for in M_PV term

" Indexed parameters - values change based on day "
p = np.genfromtxt('pv_batt_sizing_data.csv', delimiter=',', names=True)
L = p['L']        # total builing load during PV generation hours - from wx
S = p['S']        # daily load shaving requirement
M_B = p['M_B']    # dialy multiplier relate battery size to daily demand
M_PV = p['M_PV']  # daily multiplier relate PV capacity to generation - from wx

# Variables:
" Create indexes for all variables for ordering obj./constraint coeffs "
v = {}
v['x_B'] = np.array([0])
v['x_PV'] = np.array([1])
v['pu_B'] = np.arange(2,2+T)
v['po_B'] = np.arange(v['pu_B'][-1] +1,v['pu_B'][-1] +1 +T)
v['pu_PV'] = np.arange(v['po_B'][-1] +1,v['po_B'][-1] +1 +T)
v['po_PV'] = np.arange(v['pu_PV'][-1] +1,v['pu_PV'][-1] +1 +T)
v['g'] = np.arange(v['po_PV'][-1] +1,v['po_PV'][-1] +1 +T)
v['u_B'] = np.arange(v['g'][-1] +1,v['g'][-1] +1 +T)
v['e_B'] = np.arange(v['u_B'][-1] +1,v['u_B'][-1] +1 +T)
v['u_PV'] = np.arange(v['e_B'][-1] +1,v['e_B'][-1] +1 +T)
v['e_PV'] = np.arange(v['u_PV'][-1] +1,v['u_PV'][-1] +1 +T)

" Initialize vectors based on total variables length "
num_vars = sum([len(v) for v in v.values()])
c_coeff = np.zeros(num_vars)

# Objective Function:
objective = 'minimize'
c_coeff[v['e_PV']] = C_PV / (365 *Y_PV *M_PV)
c_coeff[v['u_PV']] = C_e
c_coeff[v['e_B']] = C_B / (365 *Y_B)
c_coeff[v['u_B'][0:121]] = M_B[0:121] *D_w + C_e
c_coeff[v['u_B'][121:274]] = M_B[121:274] *D_s + C_e
c_coeff[v['u_B'][274:]] = M_B[274:] *D_w + C_e

# Constraints:
constraint = {}

num_con = 1
for i in range(T):
    " Total Generation: M_PV[i] *E_PV *x_PV = g[i] "
    con_coeff = np.zeros(num_vars)
    con_coeff[v['x_PV']] = M_PV[i] *E_PV
    con_coeff[v['g'][i]] = -1
    
    constraint[num_con] = ['tot_gen_' + str(i+1), con_coeff, '=', 0]
    
    num_con = num_con +1
    
for i in range(T):
    " Battery Balance: (E_B *x_B) + u_B[i] - e_B[i] = S[i] "
    con_coeff = np.zeros(num_vars)
    con_coeff[v['x_B']] = E_B
    con_coeff[v['u_B'][i]] = 1
    con_coeff[v['e_B'][i]] = -1
    
    constraint[num_con] = ['batt_bal_' + str(i+1), con_coeff, '=', S[i]]
    
    num_con = num_con +1
    
for i in range(T):
    " PV Balance: g[i] - x_B + u_PV[i] - e_PV[i] = L[i] "
    con_coeff = np.zeros(num_vars)
    con_coeff[v['g'][i]] = 1
    con_coeff[v['x_B']] = -1
    con_coeff[v['u_PV'][i]] = 1
    con_coeff[v['e_PV'][i]] = -1
    
    constraint[num_con] = ['pv_bal_' + str(i+1), con_coeff, '=', L[i]]
    
    num_con = num_con +1 
            