# `LOCAL_PROBLEM_PLANESTRESSHARDMIXVAF_NOTAN` - solves the local problem
###  Comments
This function solves the local problem, i.e., the elastic
predictor-plastic corrector return mapping algorithm. Note that it is not
required to compute the elasto-plastic consistent tangent modulus in the
inverse problem. Assumptions: plane stress, isotropic hardening (Voce),
kinematic hardening (Armstrong-Frederick).

###  Input Arguments
`n_NR_local` (_double_) - maximum number of Newton-Raphson iterations
used for solving the local problem (return mapping algorithm)

`tol_NR_local` (_double_) - stopping tolerance for the Newton-Raphson
iteration used for solving the local problem (return mapping algorithm)

`CPlaneStress` (_double_) - elastic properties (stiffness matrix)

`SPlaneStress` (_double_) - elastic properties (compliance matrix)

`theta` (_double_) - material parameters

`H_isotropic` (_double_) - isotropic hardening parameters

`H_kinematic` (_double_) - kinematic hardening parameters

`epsilonV` (_double_) - infinitesimal strain at the current load step in
Voigt notation (epsilon_11, epsilon_22, 2*epsilon_12)

`epsilonVp_prev` (_double_) - plastic component of the infinitesimal
strain at the previous load step in Voigt notation (epsilon_p_11,
epsilon_p_22, 2*epsilon_p_12, epsilon_p_33)

`gamma_prev` (_double_) - plastic multiplier at the previous load step

`sigmaV_back_prev` (_double_) - back stress at the previous load step

###  Output Arguments
`sigmaV` (_double_) - Cauchy stress at the current load step in Voigt
notation (sigma_11, sigma_22, sigma_12)

`epsilonVp` (_double_) - plastic component of the infinitesimal
strain at the current load step in Voigt notation (epsilon_p_11,
epsilon_p_22, 2*epsilon_p_12, epsilon_p_33)

`gamma` (_double_) - plastic multiplier at the current load step

`sigma_back` (_double_) - back stress at the current load step

`converged_local` (_logical_) - indicates whether the local problem
converged

