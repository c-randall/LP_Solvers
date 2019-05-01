"""
Created on Wed Oct 31 13:17:23 2018

Author: 
    Corey R. Randall

Summary:
    Test problems used in order to check Simplex Solver coded in Python. 
    Examples include only problems with known solutions from AMPL in order to 
    check robustness of solver. Infeasible and unbounded solutions also 
    considered to ensure algorithm stops at appropriate step. Work completed
    for MEGN 586 at Colorado School of Mines during the 2018 Fall Semester.
    
"""

""" Import Useful Modules """
"-----------------------------------------------------------------------------"
import numpy as np

""" Favorite Problem - Tests a maximum problem """
"-----------------------------------------------------------------------------"
# Objective Function:
objective = 'maximize'
c_coeff = [1, 2] 

# Constraints:
constraint = {}
constraint[1] = ['name_1', [1, 1], '<=', 4]
constraint[2] = ['name_2', [1, -2], '<=', 2]
constraint[3] = ['name_3', [-2, 1], '<=', 2]

# Toggles/options for sensitivity and incoming basis selections:
sensitivity = 'on'   # runs a RHS sensitivity analysis ('on' or 'off')
incoming = 'last'    # incoming varialbe ('first' or 'last') neg c_bar value

# Tolerance Conditions: 
# Makes the (abs(x) <= tolerance) equal 0 in c_bar calculations
tolerance = 1e-10
decimals = 3

# Output:
# End Phase 1... begin Phase 2. 
#
# STOP (P2): c_bar values were all >= 0. 
#
# Solution: 
#  x = [0.667 3.333]
#
# Optimal Objective Value: 7.333

# rhs_sensitivity_report:
# {0: ('name_1_range', [2.0, 'inf'], 'dual:', 1.667),
#  1: ('name_2_range', [-6.0, 'inf'], 'dual:', -0.0),
#  2: ('name_3_range', [-6.0, 4.0], 'dual:', 0.333)}

""" Simplex Packet - Tests equations with only equalities """
"-----------------------------------------------------------------------------"
# Objective Function:
objective = 'minimize'
c_coeff = [1, -2, 0, 0, 0] 

# Constraints:
constraint = {}
constraint[1] = ['name_1', [1, 1, -1, 0, 0], '=', 2]
constraint[2] = ['name_2', [-1, 1, 0, -1, 0], '=', 1]
constraint[3] = ['name_3', [0, 1, 0, 0, 1], '=', 3]

# Toggles/options for sensitivity and incoming basis selections:
sensitivity = 'off'  # runs a RHS sensitivity analysis ('on' or 'off')
incoming = 'last'    # incoming varialbe ('first' or 'last') neg c_bar value

# Tolerance Conditions: 
# Makes the (abs(x) <= tolerance) equal 0 in c_bar calculations
tolerance = 1e-10
decimals = 3

# Output:
# End Phase 1... begin Phase 2. 
#
# STOP (P2): c_bar values were all >= 0. 
#
# Solution: 
#  x = [0. 3. 1. 2. 0.]
#
#  Optimal Objective Value: -6.0

""" Oil Company - Testing infeasibility due to constraints """
"-----------------------------------------------------------------------------"
# Objective Function:
objective = 'maximize'
c_coeff = [400, 300] 

# Constraints:
constraint = {}
constraint[1] = ['name_1', [625, 0], '<=', 6000]
constraint[2] = ['name_2', [0, 700], '<=', 5000]
constraint[3] = ['name_3', [150, 200], '>=', 12000]
constraint[4] = ['name_4', [250, 400], '>=', 20000]
constraint[5] = ['name_5', [225, 100], '>=', 15000]

# Toggles/options for sensitivity and incoming basis selections:
sensitivity = 'on'   # runs a RHS sensitivity analysis ('on' or 'off')
incoming = 'last'    # incoming varialbe ('first' or 'last') neg c_bar value

# Tolerance Conditions: 
# Makes the (abs(x) <= tolerance) equal 0 in c_bar calculations
tolerance = 1e-10
decimals = 3

# Output:
# STOP (P1): infeasible - artificials in Basis

""" Boundedness - Testing unbounded regions """
"-----------------------------------------------------------------------------"
# Objective Function:
objective = 'maximize'
c_coeff = [4, 6] 

# Constraints:
constraint = {}
constraint[1] = ['name_1', [2, -2], '<=', 6]
constraint[2] = ['name_2', [4, 0], '<=', 16]

# Toggles/options for sensitivity and incoming basis selections:
sensitivity = 'on'   # runs a RHS sensitivity analysis ('on' or 'off')
incoming = 'last'    # incoming varialbe ('first' or 'last') neg c_bar value

# Tolerance Conditions: 
# Makes the (abs(x) <= tolerance) equal 0 in c_bar calculations
tolerance = 1e-10
decimals = 3          # up to this many decimals will be reported

# Output:
# End Phase 1... begin Phase 2. 
# 
# STOP (P2): problem is unbounded

""" HW13-P2: investments with a sensitivity report """
"-----------------------------------------------------------------------------"
# Objective Function:
objective = 'maximize'
c_coeff = [13, 16, 16, 14, 39] 

# Constraints:
constraint = {}
constraint[1] = ['time_0', [11, 53, 5, 5, 29], '<=', 40]
constraint[2] = ['time_1', [3, 6, 5, 1, 34], '<=', 20]
for i in range(5):
    frac_coeff = np.zeros(5)
    frac_coeff[i] = 1
    constraint[i+3] = ['frac_'+str(i+1), frac_coeff.tolist(), '<=', 1]

# Toggles/options for sensitivity and incoming basis selections:
sensitivity = 'on'   # runs a RHS sensitivity analysis ('on' or 'off')
incoming = 'last'    # incoming varialbe ('first' or 'last') neg c_bar value

# Tolerance Conditions: 
# Makes the (abs(x) <= tolerance) equal 0 in c_bar calculations
tolerance = 1e-10
decimals = 3          # up to this many decimals will be reported

# Output:
# End Phase 1... begin Phase 2. 
#
# STOP (P2): c_bar values were all >= 0. 
# 
# Solution: 
#  x = [1.    0.201 1.    1.    0.288]
#
#  Optimal Objective Value: 57.449

# rhs_sensitivity_report:
# {0: ('time_0_range', [30.382, 78.265], 'dual:', 0.19),
#  1: ('time_1_range', [11.151, 31.276], 'dual:', 0.985),
#  2: ('frac_1_range', [0.0, 2.139], 'dual:', 7.951),
#  3: ('frac_2_range', [0.201, 'inf'], 'dual:', -0.0),
#  4: ('frac_3_range', [0.0, 2.996], 'dual:', 10.125),
#  5: ('frac_4_range', [0.0, 3.319], 'dual:', 12.063),
#  6: ('frac_5_range', [0.288, 'inf'], 'dual:', -0.0)}

""" Primal - Dual pairs: """
"-----------------------------------------------------------------------------"
"set 1..."
# Primal Objective Function:
objective = 'minimize'
c_coeff = [-1, -1] 

# Primal Constraints:
constraint = {}
constraint[1] = ['name_1', [-1, 1], '<=', 1]
constraint[2] = ['name_2', [4, -1], '<=', 10]

# Dual Objective Function:
objective = 'maximize'
c_coeff = [1, -1, 10, -10] 

# Dual Constraints:
constraint = {}
constraint[1] = ['name_1', [-1, 1, 4, -4], '<=', -1]
constraint[2] = ['name_2', [1, -1, -1, 1], '<=', -1]
constraint[3] = ['name_3', [1, -1, 0, 0], '<=', 0]
constraint[4] = ['name_4', [0, 0, 1, -1], '<=', 0]

"set 2..."
# Primal Objective Function:
objective = 'minimize'
c_coeff = [1, -2, 0, 0, 0] 

# Primal Constraints:
constraint = {}
constraint[1] = ['name_1', [1, 1, -1, 0, 0], '=', 2]
constraint[2] = ['name_2', [-1, 1, 0, -1, 0], '=', 1]
constraint[3] = ['name_3', [0, 1, 0, 0, 1], '=', 3]

# Dual Objective Function:
objective = 'maximize'
c_coeff = [2, -2, 1, -1, 3, -3] 

# Dual Constraints:
constraint = {}
constraint[1] = ['name_1', [1, -1, -1, 1, 0, 0], '<=', 1]
constraint[2] = ['name_2', [1, -1, 1, -1, 1, -1], '<=', -2]
constraint[3] = ['name_3', [-1, 1, 0, 0, 0, 0], '<=', 0]
constraint[4] = ['name_3', [0, 0, -1, 1, 0, 0], '<=', 0]
constraint[5] = ['name_3', [0, 0, 0, 0, 1, -1], '<=', 0]