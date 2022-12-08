function selection_HardMix(solutionfile)
% selection_HardMix selects a sparse solution with high fitting accuracy
%
% ## Comments
% 
% Assumptions: plane stress, isotropic hardening (Voce), kinematic
% hardening (Armstrong-Frederick).
% 
% ## Input Arguments
% 
% `solutionfile` (_char_) - filename of the solutionset
% 
% ## Output Arguments
% 
% _none_
%

addpath(genpath(pwd));
if nargin < 1
    solutionfile = '202112T01220654_RESULTS_lsqnonlin_lb_parsplit_HardMixVAF_FYS_vonMises_HardMix_plate_elliptic_holes';
end

%% Load Solution
load(solutionfile,'solution')

%% Hyperparameters
threshold_factor_cost = 1.01;
threshold_factor_parameters = 5e-3;

%% Initial Guess Solution Selection
% select from all random initial guesses the solution vector which
% corresponds to the lowest minimization function value
solution_best_random = selection_best_random_HardMixVAF(solution);

%% Lp-regularization Weighting Factor Selection
% select from all different choices of lambda_p the solution vector which
% corresponds to a high sparsity and high fitting accuracy
[parameters_select, lambda_p_select] = selection_threshold(solution_best_random.lambda_p,solution_best_random.parameters,solution_best_random.regularizer,solution_best_random.cost_without_regularizer,threshold_factor_cost,threshold_factor_parameters);

%% Print Selected Values
theta_select = parameters_select(1:end-5);
H_isotropic_select = parameters_select(end-4:end-2);
H_kinematic_select = parameters_select(end-1:end);
disp('Selected regularization factor:')
disp(['lambda_p = ', num2str(lambda_p_select)])
disp('Selected yield surface parameters:')
disp(['theta = ', num2str(theta_select)])
disp('Selected isotropic hardening parameters:')
disp(['H_isotropic = ', num2str(H_isotropic_select)])
disp('Selected kinematic hardening parameters:')
disp(['H_kinematic = ', num2str(H_kinematic_select)])

end