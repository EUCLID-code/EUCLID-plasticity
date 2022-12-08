# `LOCAL_PROBLEM_VEVPHARDMIXLIN` - solves the local problem
###  Comments
This function solves the local problem, i.e., the
viscoelastic-viscoplastic constitutive equations with mixed isotropic and
kinematic hardening.

Note:

`epsilon_vp = alpha_1 = alpha_3`

`gamma = sqrt(3/2) * alpha_2`

###  Input Arguments
`Gi` (_double_) - viscoelastic material parameters

`gi` (_double_) - viscoelastic material parameters

`Ki` (_double_) - viscoelastic material parameters

`ki` (_double_) - viscoelastic material parameters

`Ginf` (_double_) - viscoelastic material parameter

`Kinf` (_double_) - viscoelastic material parameter

`H_iso` (_double_) - isotropic hardening parameter

`H_kin` (_double_) - kinematic hardening parameter

`eta` (_double_) - viscoplastic material parameter

`sigma_0` (_double_) - yield stress

`time_inc` (_double_) - time increment

`epsilonV_2DPlaneStrain` (_double_) - infinitesimal strain at the current
load step in Voigt notation

`alphaV_prev` (_double_) - viscoelastic internal variables at the
previous load step

`epsilonVvp_prev` (_double_) - viscoplastic internal variables at the
previous load step

`gamma_prev` (_double_) - plastic multiplier at the previous load step

###  Output Arguments
`sigmaV` (_double_) - Cauchy stress at the current load step in Voigt
notation

`CC` (_double_) - consistent tangent modulus

`alphaV` (_double_) - viscoelastic internal variables at the
current load step

`epsilonVvp` (_double_) - viscoplastic internal variables at the
current load step

`gamma` (_double_) - plastic multiplier at the current load step

`viscoelastic` (_logical_) - indicates whether the load step was purely
viscoelastic

`converged_local` (_logical_) - indicates whether the local problem
converged
caution: this implementation assumes gi > 0 and ki > 0

