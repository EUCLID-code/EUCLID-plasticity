# `SELECTION_LP_THRESHOLD` - selects from all choices of `lambda_p` a solution vector with both low cost and low regularization term
###  Comments
Given multiple solutions of material parameters and their corresponding
costs and regularization terms, this function selects the solution which
has both a low cost and a low regularization term. The selection
procedure depends on a hyperparameter (`threshold_factor_cost`), which
influences how sparse the solution should be, and a hyperparameter
(`threshold_factor_parameters`), determining below which value parameters
should be discarded.

###  Input Arguments
`lambda_p` (_double_) - regularization factors

`parameters` (_double_) - multiple solutions of material parameters (every row
corresponds to one solution)

`regularizer` (_double_) - regularization terms of the given solutions
(without considering weighting factor lambda_p)

`cost_without_regularizer` (_double_) - cost function values of the given
solutions minus the regularization terms (considering weighting factor
lambda_p)

`threshold_factor_cost` (_double_) - factor for computing the cost
threshold C^{th}

`threshold_factor_parameters` (_double_) - factor for computing the
parameter threshold theta^{th}

###  Output Arguments
`parameters_select` (_double_) - selected solution of the material parameters

