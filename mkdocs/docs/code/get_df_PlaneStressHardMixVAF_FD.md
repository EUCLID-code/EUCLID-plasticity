# `GET_DF_PLANESTRESSHARDMIXVAF_FD` - evaluates the yield function and its derivatives using finite differences
###  Comments
Assumptions: plane stress, isotropic hardening (Voce), kinematic
hardening (Armstrong-Frederick).

###  Input Arguments
`theta` (_double_) - material parameters

`H_isotropic` (_double_) - isotropic hardening parameters

`sigmaV` (_double_) - Cauchy stress under plane stress conditions in
Voigt notation (sigma_11, sigma_22, sigma_22)

`gamma` (_double_) - plastic multiplier

`sigmaV_back` (_double_) - back stress

`FD` (_double_) - finite differences stepsize

`secondDerivative` (_double_) - indicates whether second derivative
should be computed

###  Output Arguments
derivatives of the yield function

