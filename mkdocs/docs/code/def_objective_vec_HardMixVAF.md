# `DEF_OBJECTIVE_VEC_HARDMIXVAF` - creates an objective function handle
###  Comments
This function creates an objective function handle. The function handle
takes the material parameters as input. Note that the function handle
does not return a scalar but a vector. The Matlab build-in optimizer
(e.g., lsqnonlin) minimizes the sum of squared entries of this vector.

A __Matlab executable file__ (.mex file) may be used to accelerate the
objective function evaluation.

###  Input Arguments
`in` (_struct_) - see `opt_input` documentation

`data` (_struct_) - FEM data

`data.results.u` (_double_) - displacement data for every load step
(every column corresponds to one load step)

`data.results.reaction` (_double_) - reaction force data for every load
step

`data.algorithm.n_NR_local` (_double_) - maximum number of Newton-Raphson
iterations used for solving the local problem (return mapping algorithm)

`data.algorithm.tol_NR_local` (_double_) - stopping tolerance for the
Newton-Raphson iteration used for solving the local problem (return
mapping algorithm)

`data.mesh.n_dof` (_double_) - number of degrees of freedom in the mesh

`data.mesh.dof` (_double_) - degrees of freedom

`data.mesh.n_element` (_double_) - number of finite elements in the mesh

`data.mesh.n_dof_per_element` (_double_) - number of degrees of freedom
per finite element

`data.mesh.element_dof` (_double_) - degrees of freedom at each finite
element element (each row provides all degrees of freedom at the
corresponding finite element)

`data.mesh.Gauss_weights` (_double_) - Gaussian quadrature weights

`data.mesh.n_Gauss_per_dim` (_double_) - number of Gauss points per
dimension

`data.mesh.detJ_GP` (_double_) - determinant of the Jacobian at each
Gauss point needed for mapping integrals between the reference and
physical finite element

`data.mesh.B_GP` (_double_) - help matrix at each Gauss point needed for
calculating shape function derivatives

`data.bc.dof_fix` (_double_) - fixed degrees of freedom

`data.bc.dof_displacement` (_double_) - free degrees of freedom

`data.bc.dof_reaction` (_double_) - degrees of freedom corresponding to
reaction forces

###  Output Arguments
`fun_vec` (*function_handle*) - objective function handle

