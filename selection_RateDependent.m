function theta_select = selection_RateDependent(filename)
% selection_RateDependent selects a sparse solution with high fitting accuracy
%
% ## Comments
% 
% _none_
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

if nargin == 0

    % noise 0
    %     filename = '202209T07202858_RESULTS_lsqnonlin_VEVP_par_perturbation_lib_VEVPHardMixLin_LE_plate_elliptic_holes.mat';
    %     filename = '202209T07203957_RESULTS_lsqnonlin_VEVP_par_perturbation_lib_VEVPHardMixLin_VE_plate_elliptic_holes.mat';
    %     filename = '202209T07204157_RESULTS_lsqnonlin_VEVP_par_perturbation_lib_VEVPHardMixLin_VEEPHardIsoLin_plate_elliptic_holes.mat';
    %     filename = '202209T07204658_RESULTS_lsqnonlin_VEVP_par_perturbation_lib_VEVPHardMixLin_LEVPHardKinLin_plate_elliptic_holes.mat';
    %     filename = '202209T07204727_RESULTS_lsqnonlin_VEVP_par_perturbation_lib_VEVPHardMixLin_VEVPHardMixLin_plate_elliptic_holes.mat';

    % noise 1
    %     filename = '202209T07205402_RESULTS_lsqnonlin_VEVP_par_perturbation_lib_VEVPHardMixLin_LE_plate_elliptic_holes_noise_1_sgolay_10.mat';
    %     filename = '202209T07205434_RESULTS_lsqnonlin_VEVP_par_perturbation_lib_VEVPHardMixLin_VE_plate_elliptic_holes_noise_1_sgolay_10.mat';
    %     filename = '202209T07210632_RESULTS_lsqnonlin_VEVP_par_perturbation_lib_VEVPHardMixLin_VEEPHardIsoLin_plate_elliptic_holes_noise_1_sgolay_10.mat';
    filename = '202209T07215606_RESULTS_lsqnonlin_VEVP_par_perturbation_lib_VEVPHardMixLin_LEVPHardKinLin_plate_elliptic_holes_noise_1_sgolay_10.mat';
    %     filename = '202209T07220636_RESULTS_lsqnonlin_VEVP_par_perturbation_lib_VEVPHardMixLin_VEVPHardMixLin_plate_elliptic_holes_noise_1_sgolay_10.mat';

    % noise 5
    %     filename = '202209T08133401_RESULTS_lsqnonlin_VEVP_par_perturbation_lib_VEVPHardMixLin_LE_plate_elliptic_holes_noise_5_sgolay_10.mat';
    %     filename = '202209T08133403_RESULTS_lsqnonlin_VEVP_par_perturbation_lib_VEVPHardMixLin_VE_plate_elliptic_holes_noise_5_sgolay_10.mat';
    %     filename = '202209T08133402_RESULTS_lsqnonlin_VEVP_par_perturbation_lib_VEVPHardMixLin_VEEPHardIsoLin_plate_elliptic_holes_noise_5_sgolay_10.mat';
    %     filename = '202209T08133402_RESULTS_lsqnonlin_VEVP_par_perturbation_lib_VEVPHardMixLin_LEVPHardKinLin_plate_elliptic_holes_noise_5_sgolay_10.mat';
    %     filename = '202209T08133403_RESULTS_lsqnonlin_VEVP_par_perturbation_lib_VEVPHardMixLin_VEVPHardMixLin_plate_elliptic_holes_noise_5_sgolay_10.mat';

end

%% Solution
load(filename,'solution')
if strcmp(solution.in.library.name,'lib_VEVPHardMixLin_v2')
    solution = compute_cost_components_VEVPHardMixLin_v2(solution);
else
    error('Not implemented.')
end

%% Data
in = opt_input_RateDependent();
load(solution.in.datafile,'data')

%% Material Library
lib_disc = feval(solution.in.library.name,solution.in.library.n_Maxwell);

%% Lowest Cost Solution
[cost_min,pos_cost_min] = min(solution.cost);

disp('True parameters:')
lib_disc.disp_theta(solution.theta_true)

%% Plot Pareto Curve
lambda_p_Pareto_idx = 1:solution.in.n_lambda_p;
theta_Pareto = solution.theta(1+end-solution.in.n_lambda_p:end,:);
lambda_p_Pareto = solution.lambda_p_set(1+end-solution.in.n_lambda_p:end);
norm_lp_Pareto = solution.norm_lp(1+end-solution.in.n_lambda_p:end);
regularizer_Pareto = solution.regularizer(1+end-solution.in.n_lambda_p:end);
cost_without_regularizer_Pareto = solution.cost_without_regularizer(1+end-solution.in.n_lambda_p:end);

plot_Pareto(lambda_p_Pareto_idx,cost_without_regularizer_Pareto,norm_lp_Pareto,...
    'labels',...
    {'$\lambda_p \cdot 10^{4}$' '$C_j$' '$\|{  \mbox{\boldmath $\theta$}  }^{\mbox{opt}}_j\|_p^p$'})
xticks([1,6,11,16,21,24])
xticklabels({'$2^0$','$2^5$','$2^{10}$','$2^{15}$','$2^{20}$','$2^{23}$'})

%% Automated Selection Algorithm
threshold_factor_cost = 1.1;
min_cost_threshold = 1e-5;
[lambda_p_select, theta_select, norm_lp_select, cost_without_regularizer_select] = ...
    selection_automatic(lambda_p_Pareto,theta_Pareto,norm_lp_Pareto,cost_without_regularizer_Pareto,threshold_factor_cost,min_cost_threshold);
lambda_p_idx_select = lambda_p_Pareto_idx(lambda_p_select == lambda_p_Pareto);

%% Threshold
theta_select = lib_disc.apply_threshold(theta_select);

disp('Automatically selected parameters:')
lib_disc.disp_theta(theta_select)

%% Plot Selection in Pareto Curve
yyaxis left
hold on
plot([lambda_p_idx_select,lambda_p_idx_select],get(gca, 'YLim'),'k--')

end



