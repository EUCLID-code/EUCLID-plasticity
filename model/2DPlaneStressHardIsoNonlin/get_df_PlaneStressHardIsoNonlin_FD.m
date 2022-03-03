function [df_dsigmaV, d2f_d2sigmaV, df_dgamma, d2f_d2gamma] = get_df_PlaneStressHardIsoNonlin_FD(theta,H_isotropic,sigmaV,gamma,FD,secondDerivative)
% get_df_PlaneStressHardIsoNonlin_FD evaluates the yield function and its
% derivatives using finite differences
%
% ## Comments
% 
% Assumptions: plane stress, isotropic hardening (Voce).
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
% `FD` (_double_) - finite differences stepsize
%
% `secondDerivative` (_double_) - indicates whether second derivative
% should be computed
%
% ## Output Arguments
% 
% derivatives of the yield function
% 

if nargin == 5
    secondDerivative = true;
end

df_dsigmaV = NaN(3,1);
d2f_d2sigmaV = NaN(3,3);
df_dgamma = NaN;
d2f_d2gamma = NaN;

% derivatives w.r.t. sigmaV
for i = 1:3
    dsigmaV = zeros(3,1);
    dsigmaV(i) = FD;
    
    % 1st derivative
    f_plus = get_f_PlaneStressHardIsoNonlin(theta,H_isotropic,sigmaV + 0.5*dsigmaV,gamma);
    f_minus = get_f_PlaneStressHardIsoNonlin(theta,H_isotropic,sigmaV - 0.5*dsigmaV,gamma);
    df_dsigmaV(i) = (f_plus-f_minus)/FD;
    % 2nd derivative
    if secondDerivative
        [df_dsigmaV_plus,~,~,~] = get_df_PlaneStressHardIsoNonlin_FD(theta,H_isotropic,sigmaV + 0.5*dsigmaV,gamma,FD,false);
        [df_dsigmaV_minus,~,~,~] = get_df_PlaneStressHardIsoNonlin_FD(theta,H_isotropic,sigmaV - 0.5*dsigmaV,gamma,FD,false);
        d2f_d2sigmaV(:,i) = (df_dsigmaV_plus - df_dsigmaV_minus)/FD;
    end
    
end

% derivatives w.r.t. gamma
% 1st derivative
f_plus = get_f_PlaneStressHardIsoNonlin(theta,H_isotropic,sigmaV,gamma + 0.5*FD);
f_minus = get_f_PlaneStressHardIsoNonlin(theta,H_isotropic,sigmaV,gamma - 0.5*FD);
df_dgamma = (f_plus-f_minus)/FD;
% 2nd derivative
if secondDerivative
    [~,~,df_dgamma_plus,~] = get_df_PlaneStressHardIsoNonlin_FD(theta,H_isotropic,sigmaV,gamma + 0.5*FD,FD,false);
    [~,~,df_dgamma_minus,~] = get_df_PlaneStressHardIsoNonlin_FD(theta,H_isotropic,sigmaV,gamma - 0.5*FD,FD,false);
    d2f_d2gamma = (df_dgamma_plus - df_dgamma_minus)/FD;
end

end