function [parameters_select, lambda_p_select] = selection_threshold(lambda_p,parameters,regularizer,cost_without_regularizer,threshold_factor_cost,threshold_factor_parameters)
% selection_Lp_threshold selects from all choices of `lambda_p` a solution vector with both low cost and low regularization term
% 
% ## Comments
% 
% Given multiple solutions of material parameters and their corresponding
% costs and regularization terms, this function selects the solution which
% has both a low cost and a low regularization term. The selection
% procedure depends on a hyperparameter (`threshold_factor_cost`), which
% influences how sparse the solution should be, and a hyperparameter
% (`threshold_factor_parameters`), determining below which value parameters
% should be discarded.
% 
% ## Input Arguments
% 
% `lambda_p` (_double_) - regularization factors
% 
% `parameters` (_double_) - multiple solutions of material parameters (every row
% corresponds to one solution)
% 
% `regularizer` (_double_) - regularization terms of the given solutions
% (without considering weighting factor lambda_p)
% 
% `cost_without_regularizer` (_double_) - cost function values of the given
% solutions minus the regularization terms (considering weighting factor
% lambda_p)
%
% `threshold_factor_cost` (_double_) - factor for computing the cost
% threshold C^{th}
%
% `threshold_factor_parameters` (_double_) - factor for computing the
% parameter threshold theta^{th}
%
% ## Output Arguments
% 
% `parameters_select` (_double_) - selected solution of the material parameters
%

cost_threshold = threshold_factor_cost*min(cost_without_regularizer); % define cost threshold above which solutions are discarded
if cost_threshold < 1e-12
    cost_threshold = 1e-12;
end
select_cost_low = find(cost_without_regularizer < cost_threshold);
lambda_p_cost_low = lambda_p(select_cost_low);
theta_cost_low = parameters(select_cost_low,:);
cost_low = cost_without_regularizer(select_cost_low);
regularizer_cost_low = regularizer(select_cost_low);

% solutions with high cost have been discarded
% find remaining solution with lowest Lp-regularizer 
regularizer_min = min(regularizer_cost_low);
select_regularizer_min = find(regularizer_cost_low == regularizer_min);
lambda_p_select = lambda_p_cost_low(select_regularizer_min);
parameters_select = theta_cost_low(select_regularizer_min,:);
if size(parameters_select,1) > 1
    % if multiple solutions have the same Lp-regularizer choose the solution with the lowest cost
    cost_regularizer_min = cost_low(select_regularizer_min);
    [~,pos] = min(cost_regularizer_min);
    lambda_p_select = lambda_p_select(pos);
    parameters_select = parameters_select(pos,:);
end

% parameter threshold
parameters_select(abs(parameters_select) < threshold_factor_parameters*parameters_select(1)) = 0;

