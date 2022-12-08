# `OPT_INPUT_RATEDEPENDENT` - provides elementary definitions needed for the optimization problem
###  Comments
Please make sure that the correct path to the FEM data is specified.

###  Input Arguments
_none_

###  Output Arguments
`in` (_struct_) - contains information needed for the optimization
problem

`in.use_mex` (_logical_) - specifies if .mex file is used for the cost
evaluation

`in.n_workers` (_double_) - number of workers used for the parallel
optimization

`in.datafile` (_char_) - filename of the FEM dataset

`in.noise_level` (_char_) - noise level to be considered (e.g., `in.noise_level=''`, `in.noise_level='1'`,
`in.noise_level='3'`, `in.noise_level='5'`)

`in.filter` (_char_) - temporal data smoothing method (e.g., `in.filter=''` for no
smoothing, `in.filter='sgolay'` for Savitzky-Golay smoothing)

`in.filter_length` (_double_) - specifies the moving window length for
the temporal smoothing

`in.library` (_obj_) - model library

`in.q_case` (_double_) - specifies which norm to use during cost function
calculation (e.g. `in.q_case=0` for L2-norm)

`in.q` (_double_) - specifies which norm to use during cost function
calculation (e.g. `in.q=2` for L2-norm)

`in.lambda_r` (_double_) - weighting factor for net reaction force
contributions

`in.lambda_p` (_double_) - weighting factor for the sparsity promoting
regularization term

`in.p` (_double_) - specifies which norm to use for the sparsity
promoting regularization

`in.lambda_p_factor` (_double_) - lambda_p is multiplied by this factor
to obtain solutions for different lambda_p

`in.n_lambda_p` (_double_) - number of different values of lambda_p to be
considered

`in.n_eval` (_double_) - evaluate cost for every n_eval-th load step

`in.method` (_char_) - method used for optimization

`in.preconditioner` (_char_) - preconditioner used for optimization

`in.preconditioner_2` (_char_) - preconditioner used for optimization

`in.batch_ratio` (_char_) - ratio defining how to split the data into two batches

`in.MaxFunEvals` (_double_) - maximum minimization function evaluations
during each optimization process

`in.TolFun` (_double_) - stopping criterion for optimization process
based on minimization function value

`in.TolX` (_double_) - stopping criterion for optimization process based
on solution values

`in.MaxTime` (_double_) - maximimum runtime of an optimization process

`in.theta_init` (_double_) - initial guess for the parameters

`in.n_guess` (_double_) - number of random initial guesses

`in.perturbation_guess_theta` (_double_) - specifies the magnitude of
perturbation applied to the preconditioning solution for generating
random initial guesses

`in.save_path` (_char_) - path for saving results

