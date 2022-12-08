# `VOLUMETRIC_DEVIATORIC_SPLITV` - applies the volumetric deviatoric split to the stress tensor in Voigt notation
###  Comments
It is `sigmaV_vol = IV_vol * sigmaV` and `sigmaV_dev = IV_dev * sigmaV`.

Assumptions:

+   Voigt notation
    
###  Input Arguments
`dim` (_double_) - dimension (default: 3)

###  Output Arguments
`IV_vol` (_double_) - matrix such that `sigmaV_vol = IV_vol * sigmaV`

`IV_dev` (_double_) - matrix such that `sigmaV_dev = IV_dev * sigmaV`

