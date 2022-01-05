[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2656688.svg)](https://doi.org/10.5281/zenodo.2590629)

## Simplex Methods
This code allows the user to solve linear programs with either the primal or dual Simplex method. These Simplex methods are algorithms that function by exploring adjacent extreme points created from the intersections of constraints. Since all linear programs are convex, the optima are located at these extreme points. Additionally, the property of convexity allows the algorithm to converge since local optima are globally optimal. 

## Objective
Software presented here is written in Python3 so that it is open and accessible to anyone who would like to use it. Problems can be formulated within the 'LP_User_Inputs.py' file by inputting the variable coefficients for the objective function and constraints. There is no need to formulate the problem in standard form since the model already incorporates a conversion routine. Constraints can be named so that the analysis of their duals is easier to interpret. It is important to note that all formulations should only use nonnegative variables. This means that any varialbes that are strictly negative should have their coefficients multiplied by -1 and the user is responsible to post process the variable value by again multiplying by -1. For variables that are unrestricted in sign, they should be formulated using x = x' - x'', where x' and x'' are strictly nonnegative. Once the final answer is found, the user can post process the output by performing the subtraction x' - x'' to get their original desired variable.

## Instructions
To install this software package onto your local machine, download all of the files from this repository. Save all of the files except for 'LP_User_Inputs.py' in the same folder and name it 'Simplex_Files' (without the quotes). The 'LP_User_Inputs.py' file should be saved in the same location as this folder (but not inside it). Make sure that your Python environment has access to numpy, scipy, sys, and os. Once these steps have been taken, the code is set up to be run and edit in any IDE with Python3 syntax. The user should not edit any routines in the 'Simplex_Files' folder. Instead, this set up has been developed so that the user should only ever need to edit the 'LP_User_Inputs.py' file in order to run any desired problems. 

Check the 'Sample_Problems.py' file for examples of how to input simple problems. For more complex problems, the user may choose to write a script that incorporates sets, parameters, variables, an objective, and constraints. An example of loading in a data file and creating an LP in a script has also been provided in the 'pv_batt_sizing.py' and 'pv_batt_sizing_sub.py' files. For more detailed information on the software package, check the 'Manual.pdf' document included with the download.

## User inputs
The 'LP_User_Inputs.py' file contains the following algorithm controls to be edited as needed by the user:
* objective = 'maximum' or 'minimum' depending on the formulation
* c_coeff = list of objective function coefficients for ALL decision variables (i.e., include zeros where necessary as placeholders)
* constraint[i] = list of constraint i characteristics as ['name', [coefficients], 'in/equality', constant]
** NOTE: The constraint indices should start with i = 1; however, as many constraints as needed can be added. Additionally, the coefficients should be in the same order as the objective coefficients listed in c_coeff (with zeros as placeholders as needed).
* sensitivity = 'on' or 'off' to either run or not run a sensitivity analysis on the right-hand-side (rhs) constants for each constraint
* pricing = 'most' or 'bland' will select the entering or exiting variable in the primal or dual simplex by selecting it based on a largest magnitude or using Bland's rule respectively
* incoming = 'first' or 'last' when using Bland's rule as the first or last indexed variable to enter or exit the basis 
* problem = 'primal' or 'dual' allows the user to solve the primal problem that they input, or to have the solver convert the problem to explicitly solve its dual instead
* method = 'primal' or 'dual' allows the user to specify whether to use the primal or dual Simplex method
* tolerance = This is an absolute tolerance for calculations in the algorithm. Any data or calculated values less than this within the magnitude of this tolerance are treated as zero.
* decimcals = The number of decimals that the output will print to the console when a converged solution is found
** NOTE: The solution is stored at a higher precision in the 'solution' dictionary output. This input simply allows for a cleaner output in the terminal when first presented with the solution.
* var_names = A list of variable names that the user wants printed at the end of the algorithm. This is only helpful for larger problems that are created with scripts. For smaller problems, set this input to 'None' (without the quotes).

## License
This tool is released under the BSD-3 clause license, see LICENSE.md for details.

## Citing the Model
 This model is versioned using Zenodo:
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2656688.svg)](https://doi.org/10.5281/zenodo.2590629)

If you use this tool as part of a scholarly work, please cite using:

> C.R. Randall (2019) Primal Simplex Method v1.2 [software]. Zenodo. https://doi.org/10.5281/zenodo.2590629

A BibTeX entry for LaTeX users is

```TeX
@misc{LP_Solver,
    author = {Corey R. Randall},
    year = 2019,
    title = {Python Simplex Solver v1.2},
    doi = {10.5281/zenodo.2590629},
    url = {https://github.com/c-randall/LP_Solver},
}
```

In both cases, please update the entry with the version used. The DOI for the latest version is
given in the badge at the top, or alternately <https://doi.org/10.5281/zenodo.2590629> will
take you to the latest version (and generally represents all versions).
