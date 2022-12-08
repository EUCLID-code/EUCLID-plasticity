# `SELECTION_AUTOMATIC` - selects from all choices of `lambda_p` a solution vector with both low cost and low lp-(semi)-norm
###  Comments
Given multiple solutions of material parameters and their corresponding
costs and lp-(semi)-norm values, this function selects the solution which
has both a low cost and a low lp-(semi)-norm. The selection
procedure depends on a hyperparameter (`threshold_factor_cost`), which
influences how sparse the solution should be, and a hyperparameter
(`threshold_factor_parameters`), determining below which value parameters
should be discarded.

###  Input Arguments
`lambda_p` (_double_) - regularization factors

`parameters` (_double_) - multiple solutions of material parameters (every row
corresponds to one solution)

`norm_lp` (_double_) - lp-(semi)-norm values (without considering
weighting factor lambda_p)

`cost_without_regularizer` (_double_) - cost function values of the given
solutions minus the regularization terms (considering weighting factor
lambda_p)

`threshold_factor_cost` (_double_) - factor for computing the cost
threshold C^{th}

`min_cost_threshold` (_double_) - minimum value of the cost threshold
C^{th}

###  Output Arguments
`lambda_p_select` (_double_) - selected regularization factor

`parameters_select` (_double_) - selected solution of the material parameters

`norm_lp_select` (_double_) - selected lp-(semi)-norm value

`cost_without_regularizer_select` (_double_) - selected cost function value

