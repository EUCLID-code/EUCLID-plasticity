function f = get_f_PlaneStressHardMixVAF(theta,H_isotropic,sigmaV,gamma,sigmaV_back)
% get_f_PlaneStressHardMixVAF evaluates the yield function
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
% `f` (_double_) - value of the yield function

sigmaV_rel = sigmaV - sigmaV_back; % relative stress = stress - back stress

% the yield function under kinematic hardening is calculated by evaluating
% the yield function for the relative stress
f = get_f_PlaneStressHardIsoNonlin(theta,H_isotropic,sigmaV_rel,gamma); 

end