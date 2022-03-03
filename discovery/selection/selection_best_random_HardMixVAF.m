function solution_best_random = selection_best_random_HardMixVAF(solution)
% selection_best_random_HardIsoNonlin selects from a set of solutions for each choice of `lambda_p` the lowest cost solution
%
% ## Comments
% 
% Note that `cost` denotes here the total minimization function, i.e.,
% including the regularization term. `cost_without_regularizer` denotes the
% cost without the weighted regularization term.
% 
% ## Input Arguments
% 
% `solution` (_struct_) - set of solutions
% 
% ## Output Arguments
% 
% `solution_best_random` (_struct_) - lowest cost solution for each choice
% of `lambda_p`

solution.parameters = [solution.theta, solution.H_isotropic, solution.H_kinematic];
lambda_p_values = unique(solution.lambda_p_set); % different values of lambda_p
n_different_lambda_p = length(lambda_p_values) ; % number of different values of lambda_p

parameters_best_random = zeros(n_different_lambda_p,size(solution.parameters,2)); % allocation
cost_best_random = zeros(n_different_lambda_p,1); % allocation
regularizer_best_random = zeros(n_different_lambda_p,1); % allocation
cost_without_regularizer_best_random = zeros(n_different_lambda_p,1); % allocation

for idx_lambda_p = 1:n_different_lambda_p % loop over different values of lambda_p
    % consider only solutions that correspond to a certain value of lambda_p
    consider = (solution.lambda_p_set == lambda_p_values(idx_lambda_p));
    parameters_consider = solution.parameters(consider,:);
    cost_consider = solution.cost(consider);
    % find lowest cost solution corresponding to a certain value of lambda_p
    [~,pos] = min(cost_consider);
    parameters_best_random(idx_lambda_p,:) = parameters_consider(pos,:);
    cost_best_random(idx_lambda_p) = cost_consider(pos);
    % compute the regularization term for the lowest cost solution
    regularizer_best_random(idx_lambda_p) = sum(abs(parameters_best_random(idx_lambda_p,2:end-5)).^solution.in.p);
    % compute the cost without the weighted regularization term (lambda_p*||theta||_p^p)
    cost_without_regularizer_best_random(idx_lambda_p) = cost_best_random(idx_lambda_p) - lambda_p_values(idx_lambda_p)*regularizer_best_random(idx_lambda_p);
end

solution_best_random.lambda_p = lambda_p_values;
solution_best_random.parameters = parameters_best_random;
solution_best_random.cost = cost_best_random;
solution_best_random.regularizer = regularizer_best_random;
solution_best_random.cost_without_regularizer = cost_without_regularizer_best_random;






