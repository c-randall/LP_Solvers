[![DOI](https://zenodo.org/badge/160264638.svg)](https://zenodo.org/badge/latestdoi/160264638)

## Simplex Methods
This code allows the user to solve linear programs with either the primal or dual Simplex method. Theses Simplex methods are algorithms that function by exploring adjacent extreme points created from the intersections of constraints. Since all linear programs are convex, the optima are located at these extreme points. Additionally, the property of convexity allows the algorithm to converge since local optima are globally optimal. 

## Objective
Software presented here is written in Python3 so that it is open and accessible to anyone who would like to use it. Problems can be formulated within the LP_User_Inputs.py file by inputting the variable coefficients for the objective function and constraints. There is no need to formulate the problem in standard form since the model already incorporates a conversion file. Constraints can be named so that the analysis of their duals is easier to interpret. 

## Instructions
Use the Sample_Problems.py file for examples of how to input simple problems. For more complex problems, the user may choose to write a script that incorporates sets, parameters, varialbes, an objective, and constraints. An example of loading in a data file and creating an LP in a script has also been provided in the pv_batt_sizing.py and pv_batt_sizing_sub.py files.

## User inputs
The LP_User_Inputs.py file contains the following algorithm controls to be edited as needed by the user:
* objective = 'maximum' or 'minimum' depending on the formulation
* c_coeff = list of objective function coefficients for ALL decision variables (i.e., include zeros where necessary as placeholders)
* constraint[i] = list of constraint i characteristics as ['name', [coefficients], 'in/equality', constant]
** NOTE: The constraint indices should start with i = 1; however, as many constraints as needed can be added. Additionally, the coefficients should be in the same order as the objective coefficients listed in c_coeff (with zeros as placeholders as needed).
* sensitivity = 'on' or 'off' to either run or not run a sensitivity analysis on the right-hand-side (rhs) constants for each constraint
* pricing = 'most' or 'bland' will select the entering or exiting variable in the primal or dual simplex by selecting it based on a largest magnitude or using Bland's rule respectively
* incoming = 'first' or 'last' when using Bland's rule as the first or last indexed variable to enter the basis 
* problem = 'primal' or 'dual' allows the user to solve the primal problem that they input, or to have the solver convert the problem to explicitly solve the dual instead
* method = 'primal' or 'dual' allows the user to specify whether to use the primal or dual Simplex method
* tolerance = This is an absolute tolerance for calculations in the algorithm. Any calculated values less than this tolerance are treated as zero.
* decimcals = The number of decimals that the output will print out when a converged solution is found
** NOTE: The solution is stored at a higher precision in the 'solution' dictionary output. This input simply allows for a cleaner output in the terminal when first presented with the solution.
* var_names = A list of variable names that the user wants printed at the end of the algorithm. This is only helpful for larger problems that are created with scripts. For smaller problems, this set this to an empty list, i.e. [].

## License
This tool is released under the BSD-3 clause license, see LICENSE.md for details.

## Citing the Model
 This model is versioned using Zenodo:
[![DOI](https://zenodo.org/badge/160264638.svg)](https://zenodo.org/badge/latestdoi/160264638)

If you use this tool as part of a scholarly work, please cite using:

> C.R. Randall (2019) Primal Simplex Method v1.0 [software]. Zenodo. https://doi.org/10.5281/zenodo.2590630

A BibTeX entry for LaTeX users is

```TeX
@misc{2dPorousFlux,
    author = {Corey R. Randall},
    year = 2019,
    title = {Primal Simplex Method v1.0},
    doi = {10.5281/zenodo.1317600},
    url = {https://github.com/c-randall/Simplex-Method},
}
```

In both cases, please update the entry with the version used. The DOI for the latest version is
given in the badge at the top, or alternately <https://doi.org/10.5281/zenodo.2590630> will
take you to the latest version (and generally represents all versions).
