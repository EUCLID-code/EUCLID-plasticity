function solution = compute_cost_components_VEVPHardMixLin_v2(solution)

% compute the lp-(semi)-norm
solution.norm_lp = zeros(size(solution.cost));
for idx = 1:length(solution.norm_lp)
    solution.norm_lp(idx) = sum(abs(solution.theta(idx,3:end)).^solution.in.p);
end
% compute the regularization term
solution.regularizer = solution.lambda_p_set.*solution.norm_lp;
% compute the cost excluding the regularization term
solution.cost_without_regularizer = solution.cost - solution.regularizer;

end
