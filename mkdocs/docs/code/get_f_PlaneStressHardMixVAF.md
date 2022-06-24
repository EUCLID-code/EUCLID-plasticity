# `GET_F_PLANESTRESSHARDMIXVAF` - evaluates the yield function
###  Comments
Assumptions: plane stress, isotropic hardening (Voce), kinematic
hardening (Armstrong-Frederick).

###  Input Arguments
`theta` (_double_) - material parameters

`H_isotropic` (_double_) - isotropic hardening parameters

`sigmaV` (_double_) - Cauchy stress under plane stress conditions in
Voigt notation (sigma_11, sigma_22, sigma_12)

`gamma` (_double_) - plastic multiplier

`sigmaV_back` (_double_) - back stress

###  Output Arguments
`f` (_double_) - value of the yield function

