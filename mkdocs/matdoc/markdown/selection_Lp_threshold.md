# `SELECTION_LP_THRESHOLD` - selects from all choices of `lambda_p` a solution vector with both low cost and low regularization term
###  Comments
Given multiple solutions of material parameters and their corresponding
costs and regularization terms, this function selects the solution which
has both a low cost and a low regularization term. The selection
procedure depends on one hyperparameter (`threshold_factor`), which
influences how sparse the solution should be.

###  Input Arguments
`theta` (_double_) - multiple solutions of material parameters (every row
corresponds to one solution)

`cost` (_double_) - cost function values of the given solutions

`regularizer` (_double_) - regularization term of the given solutions
(without considering weighting factor lambda_p)

`threshold_factor` (_double_) - factor for computing the cost threshold
C^{th}

###  Output Arguments
`theta_select` (_double_) - selected solution of the material parameters

