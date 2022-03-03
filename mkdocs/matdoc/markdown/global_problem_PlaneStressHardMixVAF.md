# `GLOBAL_PROBLEM_PLANESTRESSHARDMIXVAF` - solves the global finite element problem
###  Comments
This function solves the global finite element problem, i.e., the
residuals of the discretized weak formulation are minimized. Assumptions:
plane stress, isotropic hardening (Voce), kinematic hardening
(Armstrong-Frederick).

###  Input Arguments
`n_save_YieldSurface` (_double_) - 

`n_NR_global` (_double_) - maximum number of Newton-Raphson iterations
used for solving the global problem

`tol_NR_global` (_double_) - stopping tolerance for the Newton-Raphson
iteration used for solving the global problem

`n_NR_local` (_double_) - maximum number of Newton-Raphson iterations
used for solving the local problem (return mapping algorithm)

`tol_NR_local` (_double_) - stopping tolerance for the Newton-Raphson
iteration used for solving the local problem (return mapping algorithm)

`CPlaneStress` (_double_) - elastic properties (stiffness matrix)

`SPlaneStress` (_double_) - elastic properties (compliance matrix)

`theta` (_double_) - material parameters

`H_isotropic` (_double_) - isotropic hardening parameters

`H_kinematic` (_double_) - kinematic hardening parameters

`n_dof` (_double_) - see `mesh_info`

`n_element` (_double_) - see `mesh_info`

`n_dof_per_element` (_double_) - see `mesh_info`

`element_dof` (_double_) - see `mesh_info`

`Gauss_weights` (_double_) - see `mesh_info`

`n_Gauss_per_dim` (_double_) - see `mesh_info`

`detJ_GP` (_double_) - see `mesh_info`

`B_GP` (_double_) - see `mesh_info`

`dof_fix` (_double_) - fixed degrees of freedom

`dof_displacement` (_double_) - degrees of freedom under displacement
control

`dof_reaction` (_double_) - degrees of freedom associated with the
reaction forces

`max_displacement` (_double_) - maximum applied displacement

`n_step` (_double_) - number of load steps

`displacement_increment` (_double_) - applied displacement increment

###  Output Arguments
`u` (_double_) - displacement at each node for each load step

`reaction` (_double_) - reaction forces for each load step

`epsilonVp_GP` (_double_) - plastic component of the infinitesimal
strain at each Gauss point for the last load step

`gamma_GP` (_double_) - plastic multiplier at each Gauss point for the last load step

`sigmaV_back_GP` (_double_) - back stress at each Gauss point for the last load step

`YieldSurface` (_double_) - stresses at which yielding has occured 

