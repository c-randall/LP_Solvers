[![DOI](https://zenodo.org/badge/160264638.svg)](https://zenodo.org/badge/latestdoi/160264638)

# Primal Simplex Method
The primal Simplex method is an algorithm that can be used to solve linear optimization problems. It functions by exploring adjacent extreme points created from the intersections of constraints. Since all linear programs are convex, the optima are located at these extreme points. Additionally, the property of convexity allows the algorithm to converge since local optima are globally optimal. 

## Objective
Software presented here is written in Python3 so that it is open and accessible to anyone who would like to use it. Problems can be formulated within the LP_runner.py file by inputting the variable coefficients for the objective function and constraints. There is no need to formulate the problem in standard form since the model already incorporates a conversion file. Constraints can be named so that the analysis of their duals is easier to interpret. Use the Sample_Problems.py file for examples of how to input problems into the runner.

## User inputs
The LP_Runner.py file contains the following variables to be edited as needed by the user:
* objective = 'maximum' or 'minimum' depending on the formulation
* c_coeff = list of objective function coefficients for ALL decision variables (i.e., include zeros where necessary as placeholders)
* constraint[i] = list of constraint i characteristics as ['name', [coefficients], 'in/equality', constant]
** NOTE: The constraint indices should start with i = 1; however, as many constraints as needed can be added. Additionally, the coefficients should be in the same order as the objective coefficients listed in c_coeff (with zeros as placeholders as needed).
* sensitivity = 'on' or 'off' to either run or not run a sensitivity analysis on the right-hand-side (rhs) constants for each constraint
* incoming = 'first' or 'last' to use Bland's rule as the first or last indexed variable to enter the basis 
** NOTE: Changing this can help with time and stability issues for larger problems since a more complex pricing scheme has not yet been written.
* tolerance = small value such that any reduced cost within this tolerance is treated as zero due to issues with finite computational precision 
* decimcals = number of decimals that the output will print out when a converged solution is found
** NOTE: The solution is stored at a higher precision in the 'solution' dictionary output. This input simply allows for a cleaner output in the terminal when first presented with the solution.

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
