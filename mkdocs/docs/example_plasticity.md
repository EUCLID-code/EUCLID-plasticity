# Example
It follows a step-by-step description of the implemented algorithm for an exemplary dataset.
In particular, EUCLID (Efficient Unsupervised Constitutive Law Identification & Discovery) is applied to the displacement and net reaction force data
that were generated based on the material model F1.
The goal is to use the data to discover the material model, without knowing its specific mathematical form a priori.

## Data and Parameters
The input data for EUCLID and the parameters for the optimization process are defined in `opt_input.m`.
First, it is defined whether a __Matlab executable file__ (.mex file) should be used to accelerate the objective function evaluation (see [__Matlab Coder__](https://mathworks.com/help/coder/) and [__MEX File Functions__](https://mathworks.com/help/matlab/call-mex-file-functions.html))
and over how many __Matlab workers__ the independent optimization processes should be distributed (see [__Parallel Computing Toolbox__](https://mathworks.com/products/parallel-computing.html) and [__parfor__](https://mathworks.com/help/parallel-computing/parfor.html)).
```Matlab
in.use_mex = true;
in.n_workers = 24;
```
To consider the data corresponding to material model F1, the variable `dataID` has to be set to three.
```Matlab
dataID = 3;
```
To also consider noise and temporal denoising, we set:
```Matlab
in.noise_level = '1';
in.filter = 'sgolay';
in.filter_length = 50;
```
According to these choices the following string is assigned to the variable `in.datafile`:
```Matlab
'20210518T174453_DATA_FYS_1_plate_elliptic_holes_noise_1_sgolay_50'
```
The string contains a timestamp that indicates when the data were generated and contains information about the material model, the geometry, the noise level and the smoothing method, respectively.
To load the .mat datafile with the above name, a path has to be added by changing the variable `pathname` in `opt_input.m` accordingly.
The .mat datafile contains a _struct_ array called `data` with all information provided by the finite element simulation,
such as for example the nodal displacements in `data.results.u` (see animation below) and net reaction forces in `data.results.reaction`.

<img src="/img/deformation_u_F1_noise_1_sgolay_50.gif" alt="Deformation" width="500"/>

<sub>Animation 1: Magnitude of displacement obtained from finite element simulations. The displacement is magnified by a factor of 100.</sub>

In `opt_input.m` we further define parameters and hyperparameters for the optimization problem.
All input information needed for the optimization problem are stored in a _struct_ array called `in`.

## Optimization
After defining the elementary inputs in `opt_input.m`, we start the optimization procedure by running `opt.m`.
In `opt.m` the parameters and the data are loaded and (as mentioned above) stored in the _struct_ arrays `in` and `data`, respectively.
```Matlab
in = opt_input();
load(in.datafile,'data')
```

<img src="/img/Fourier_flexibility.gif" alt="Yield Surface" width="400"/>

<sub>Animation 2: Yield surface for different choices of `theta`.</sub>

The goal is now to find the unknown material parameters `theta`, which influence the shape of the plastic yield surface (see animation above) and consequently govern the material behavior.
To find the material parameters, an objective *function_handle* `fun_vec` is defined. 
```Matlab
fun_vec = def_objective_vec(in,data);
```
The *function_handle* takes the material parameters as input and returns a vector, whose sum of squared entries equals the objective function.
Note that the *function_handle* returns a vector instead of a scalar as this is a requirement for using the Matlab build in optimizer `lsqnonlin` (see below). 

Taking a closer look at `objective_vec.m` provides a deeper insight in the internal processes of the objective function calculation.
Here, nested loops over all load steps, finite elements and quadrature points are implemented for calculating the internal forces.
To this end, the elasto-plastic return mapping algorithm, whose behavior depends on the yield function and hence the material parameters `theta`, needs to be applied (see `local_problem_PlaneStress_notan.m`).
For decreasing the computational time needed for the cost function evalution, the function `objective_vec.m` is converted into a __Matlab executable file__ (.mex file).

After defining the objective function we apply a preconditioning optimization by minimizing the objective function with respect to the first material parameter while keeping the other parameters fixed.
```Matlab
[theta_init(1), min_cost, ~, ~, ~] = lsqnonlin(fun_vec,theta_init(1),[],[],options);
```
In this way, the approximate size of the yield surface is estimated.
Based on the preconditioning solution, multiple random initial guesses are generated and the objective function is minimized for different choices of the hyperparameter `lambda_p`.
The solutions and their corresponding objective function values are stored in a _struct_ array called `solution`, which is saved as a .mat file with the name:
```Matlab
`202106T29085533_RESULTS_lsqnonlin_Pareto_random_par_FYS_1_plate_elliptic_holes_noise_1_sgolay_50`
```

## Model Selection
As a final postprocessing step, a model with both high parsimony and low cost function value needs to be selected from the generated set of solutions.
To this end, `selection.m` is executed.
Here, the hyperparameters for the model selection are chosen as:
```Matlab
threshold_factor_C = 1.01;
threshold_factor_theta = 5e-3;
```
As solutions were generated based on multiple initial guesses,
we choose for every choice of `lambda_p` the initial guess that corresponds to the lowest objective function value, resulting in one solution for every choice of `lambda_p`.
Afterwards, solutions with cost function values above a threshold are discarded as they are expected to provide low fitting accuracies.
From the remaining solutions, the solution with the highest sparsity, i.e., with the lowest regularization term, is selected.
Lastly, material parameters below a threshold are discarded.

The final solution is: 
```Matlab
>> selection('202106T29085533_RESULTS_lsqnonlin_Pareto_random_par_FYS_1_plate_elliptic_holes_noise_1_sgolay_50')
Selected Solution:
    0.2187    0.0220         0         0         0
```
The algorithm correctly identifies which features have to be active to describe the material behavior and at the same time correctly removes irrelevant features.
The correct mathematical form of the material model F1 is identified,
however, the material parameter values vary slightly compared to the true parameters as a consequence of the noise added to the displacement data.

