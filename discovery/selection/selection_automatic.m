function [lambda_p_select, parameters_select, norm_lp_select, cost_without_regularizer_select] = selection_automatic(lambda_p,parameters,norm_lp,cost_without_regularizer,threshold_factor_cost,min_cost_threshold)
% selection_automatic selects from all choices of `lambda_p` a solution vector with both low cost and low lp-(semi)-norm
%
% ## Comments
%
% Given multiple solutions of material parameters and their corresponding
% costs and lp-(semi)-norm values, this function selects the solution which
% has both a low cost and a low lp-(semi)-norm. The selection
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
% `norm_lp` (_double_) - lp-(semi)-norm values (without considering
% weighting factor lambda_p)
%
% `cost_without_regularizer` (_double_) - cost function values of the given
% solutions minus the regularization terms (considering weighting factor
% lambda_p)
%
% `threshold_factor_cost` (_double_) - factor for computing the cost
% threshold C^{th}
%
% `min_cost_threshold` (_double_) - minimum value of the cost threshold
% C^{th}
%
% ## Output Arguments
%
% `lambda_p_select` (_double_) - selected regularization factor
%
% `parameters_select` (_double_) - selected solution of the material parameters
%
% `norm_lp_select` (_double_) - selected lp-(semi)-norm value
%
% `cost_without_regularizer_select` (_double_) - selected cost function value
%

cost_threshold = threshold_factor_cost*min(cost_without_regularizer); % define cost threshold above which solutions are discarded
if cost_threshold < min_cost_threshold
    cost_threshold = min_cost_threshold;
end
select_cost_low = find(cost_without_regularizer < cost_threshold);
lambda_p_cost_low = lambda_p(select_cost_low);
parameters_cost_low = parameters(select_cost_low,:);
norm_lp_cost_low = norm_lp(select_cost_low);
cost_without_regularizer_cost_low = cost_without_regularizer(select_cost_low);

% solutions with high cost have been discarded
% find remaining solution with lowest lp-(semi)-norm value
norm_lp_min = min(norm_lp_cost_low);
select_regularizer_min = find(norm_lp_cost_low == norm_lp_min);
lambda_p_select = lambda_p_cost_low(select_regularizer_min);
parameters_select = parameters_cost_low(select_regularizer_min,:);
norm_lp_select = norm_lp_cost_low(select_regularizer_min);
cost_without_regularizer_select = cost_without_regularizer_cost_low(select_regularizer_min);
if size(parameters_select,1) > 1
    % if multiple solutions have the same lp-(semi)-norm value choose the solution with the lowest cost
    cost_regularizer_min = cost_without_regularizer_cost_low(select_regularizer_min);
    [~,pos] = min(cost_regularizer_min);
    lambda_p_select = lambda_p_select(pos);
    parameters_select = parameters_select(pos,:);
    norm_lp_select = norm_lp_select(pos,:);
    cost_without_regularizer_select = cost_without_regularizer_select(pos,:);
end

