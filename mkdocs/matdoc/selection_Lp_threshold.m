function theta_select = selection_Lp_threshold(theta,cost,regularizer,threshold_factor)
% selection_Lp_threshold selects from all choices of `lambda_p` a solution vector with both low cost and low regularization term
%
% ## Comments
% 
% Given multiple solutions of material parameters and their corresponding
% costs and regularization terms, this function selects the solution which
% has both a low cost and a low regularization term. The selection
% procedure depends on one hyperparameter (`threshold_factor`), which
% influences how sparse the solution should be.
% 
% ## Input Arguments
% 
% `theta` (_double_) - multiple solutions of material parameters (every row
% corresponds to one solution)
% 
% `cost` (_double_) - cost function values of the given solutions
%
% `regularizer` (_double_) - regularization term of the given solutions
% (without considering weighting factor lambda_p)
% 
% `threshold_factor` (_double_) - factor for computing the cost threshold
% C^{th}
%
% ## Output Arguments
% 
% `theta_select` (_double_) - selected solution of the material parameters
%

cost_threshold = threshold_factor*min(cost); % define cost threshold above which solutions are discarded
select_cost_low = find(cost < cost_threshold);
theta_cost_low = theta(select_cost_low,:);
cost_low = cost(select_cost_low);
regularizer_cost_low = regularizer(select_cost_low);

% solutions with high cost have been discarded
% find remaining solution with lowest Lp-regularizer 
regularizer_min = min(regularizer_cost_low);
select_regularizer_min = find(regularizer_cost_low == regularizer_min);
theta_select = theta_cost_low(select_regularizer_min,:);
if size(theta_select,1) > 1
    % if multiple solutions have the same Lp-regularizer choose the solution with the lowest cost
    cost_regularizer_min = cost_low(select_regularizer_min);
    [~,pos] = min(cost_regularizer_min);
    theta_select = theta_select(pos,:);
end

end