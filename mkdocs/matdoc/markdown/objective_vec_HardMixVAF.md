# `OBJECTIVE_VEC_HARDMIXVAF` - calculates the objective minimization function
###  Comments 
This function defines the objective function that depends on
the given data and on theta. Note that the function
does not return a scalar but a vector. The Matlab build-in optimizer
(e.g., lsqnonlin) minimizes the sum of squared entries of this vector.

This function may be converted into a __Matlab executable file__ (.mex
file) with the [__Matlab Coder__](https://mathworks.com/help/coder/) to
accelerate the objective function evaluation.

###  Input Arguments
`paramters` (_double_) - material parameters

`q_case` (_double_) - see `opt_input` documentation

`q` (_double_) - see `opt_input` documentation

`lambda_r` (_double_) - see `opt_input` documentation

`lambda_p` (_double_) - see `opt_input` documentation

`p` (_double_) - see `opt_input` documentation

`n_eval` (_double_) - see `opt_input` documentation

`u` (_double_) - displacement data for every load step (every column
corresponds to one load step)

`reaction` (_double_) - reaction force data for every load step

`n_NR_local` (_double_) - maximum number of Newton-Raphson iterations
used for solving the local problem (return mapping algorithm)

`tol_NR_local` (_double_) - stopping tolerance for the Newton-Raphson
iteration used for solving the local problem (return mapping algorithm)

`CPlaneStress` (_double_) - known elastic properties (stiffness matrix)

`SPlaneStress` (_double_) - known elastic properties (compliance matrix)

`n_dof` (_double_) - number of degrees of freedom in the mesh

`n_element` (_double_) - number of finite elements in the mesh

`n_dof_per_element` (_double_) - number of degrees of freedom per finite
element

`element_dof` (_double_) - degrees of freedom at each finite element
element (each row provides all degrees of freedom at the corresponding
finite element)

`Gauss_weights` (_double_) - Gaussian quadrature weights

`n_Gauss_per_dim` (_double_) - number of Gauss points per dimension

`detJ_GP` (_double_) - determinant of the Jacobian at each Gauss point
needed for mapping integrals between the reference and physical finite
element

`B_GP` (_double_) - help matrix at each Gauss point needed for
calculating shape function derivatives

`dof_free` (_double_) - free degrees of freedom

`dof_reaction` (_double_) - degrees of freedom corresponding to reaction
forces

###  Output Arguments
`cost` (_double_) - vector, whose sum of squared entries must be
minimized (note that the regularization term is included as the last
entry of this vector)

