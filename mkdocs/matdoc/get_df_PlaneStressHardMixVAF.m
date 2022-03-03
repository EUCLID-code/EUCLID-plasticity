function [f, df_dsigmaV, d2f_d2sigmaV, df_dgamma, d2f_d2gamma, d2f_dsigmaVdgamma, df_dsigmaV_back, d2f_d2sigmaV_back, d2f_dsigmaVdsigmaV_back] = get_df_PlaneStressHardMixVAF(theta,H_isotropic,sigmaV,gamma,sigmaV_back)
% get_df_PlaneStressHardMixVAF evaluates the yield function and its derivatives
%
% ## Comments
% 
% Assumptions: plane stress, isotropic hardening (Voce), kinematic
% hardening (Armstrong-Frederick).
% 
% ## Input Arguments
% 
% `theta` (_double_) - material parameters
% 
% `H_isotropic` (_double_) - isotropic hardening parameters
%
% `sigmaV` (_double_) - Cauchy stress under plane stress conditions in
% Voigt notation (sigma_11, sigma_22, sigma_22)
% 
% `gamma` (_double_) - plastic multiplier
%
% `sigmaV_back` (_double_) - back stress
%
% ## Output Arguments
% 
% yield function value and derivatives of the yield function
% 

sigmaV_rel = sigmaV - sigmaV_back; % relative stress = stress - back stress

% the yield function under kinematic hardening is calculated by evaluating
% the yield function for the relative stress
% [f, df_dsigmaV_rel, d2f_d2sigmaV_rel] = get_df_PlaneStress(theta,sigmaV_rel);
[f, df_dsigmaV, d2f_d2sigmaV, df_dgamma, d2f_d2gamma, d2f_dsigmaVdgamma] = get_df_PlaneStressHardIsoNonlin(theta,H_isotropic,sigmaV_rel,gamma);

df_dsigmaV_back = - df_dsigmaV;
d2f_d2sigmaV_back = d2f_d2sigmaV;
d2f_dsigmaVdsigmaV_back = - d2f_d2sigmaV;

end