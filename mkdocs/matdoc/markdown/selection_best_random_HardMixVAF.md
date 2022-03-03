# `SELECTION_BEST_RANDOM_HARDISONONLIN` - selects from a set of solutions for each choice of `lambda_p` the lowest cost solution
###  Comments
Note that `cost` denotes here the total minimization function, i.e.,
including the regularization term. `cost_without_regularizer` denotes the
cost without the weighted regularization term.

###  Input Arguments
`solution` (_struct_) - set of solutions

###  Output Arguments
`solution_best_random` (_struct_) - lowest cost solution for each choice
of `lambda_p`

