function [df_dsigmaV, d2f_d2sigmaV, df_dgamma, d2f_d2gamma, d2f_dsigmaVdgamma, df_dsigmaV_back, d2f_d2sigmaV_back, d2f_dsigmaVdsigmaV_back] = get_df_PlaneStressHardMixVAF_FD(theta,H_isotropic,sigmaV,gamma,sigmaV_back,FD,secondDerivative)
% get_df_PlaneStressHardMixVAF_FD evaluates the yield function and its derivatives using finite differences
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
% `FD` (_double_) - finite differences stepsize
%
% `secondDerivative` (_double_) - indicates whether second derivative
% should be computed
% 
% ## Output Arguments
% 
% derivatives of the yield function
% 

if nargin == 6
    secondDerivative = true;
end

df_dsigmaV = NaN(3,1);
d2f_d2sigmaV = NaN(3,3);
df_dgamma = NaN;
d2f_d2gamma = NaN;
d2f_dsigmaVdgamma = NaN(3,1);
df_dsigmaV_back = NaN(3,1);
d2f_d2sigmaV_back = NaN(3,3);
d2f_dsigmaVdsigmaV_back = NaN(3,3);

% derivatives w.r.t. sigmaV
for i = 1:3
    dsigmaV = zeros(3,1);
    dsigmaV(i) = FD;
    
    % 1st derivative
    f_plus = get_f_PlaneStressHardMixVAF(theta,H_isotropic,sigmaV + 0.5*dsigmaV,gamma,sigmaV_back);
    f_minus = get_f_PlaneStressHardMixVAF(theta,H_isotropic,sigmaV - 0.5*dsigmaV,gamma,sigmaV_back);
    df_dsigmaV(i) = (f_plus-f_minus)/FD;
    % 2nd derivative
    if secondDerivative
        df_dsigmaV_plus = get_df_PlaneStressHardMixVAF_FD(theta,H_isotropic,sigmaV + 0.5*dsigmaV,gamma,sigmaV_back,FD,false);
        df_dsigmaV_minus = get_df_PlaneStressHardMixVAF_FD(theta,H_isotropic,sigmaV - 0.5*dsigmaV,gamma,sigmaV_back,FD,false);
        d2f_d2sigmaV(:,i) = (df_dsigmaV_plus - df_dsigmaV_minus)/FD;
    end
    
end

% derivatives w.r.t. gamma
% 1st derivative
f_plus = get_f_PlaneStressHardMixVAF(theta,H_isotropic,sigmaV,gamma + 0.5*FD,sigmaV_back);
f_minus = get_f_PlaneStressHardMixVAF(theta,H_isotropic,sigmaV,gamma - 0.5*FD,sigmaV_back);
df_dgamma = (f_plus-f_minus)/FD;
% 2nd derivative
if secondDerivative
    [~,~,df_dgamma_plus] = get_df_PlaneStressHardMixVAF_FD(theta,H_isotropic,sigmaV,gamma + 0.5*FD,sigmaV_back,FD,false);
    [~,~,df_dgamma_minus] = get_df_PlaneStressHardMixVAF_FD(theta,H_isotropic,sigmaV,gamma - 0.5*FD,sigmaV_back,FD,false);
    d2f_d2gamma = (df_dgamma_plus - df_dgamma_minus)/FD;
end

% derivatives w.r.t. sigmaV_back and mixed derivative
for i = 1:3
    dsigmaV_back = zeros(3,1);
    dsigmaV_back(i) = FD;
    
    % 1st derivative
    f_plus = get_f_PlaneStressHardMixVAF(theta,H_isotropic,sigmaV,gamma,sigmaV_back + 0.5*dsigmaV_back);
    f_minus = get_f_PlaneStressHardMixVAF(theta,H_isotropic,sigmaV,gamma,sigmaV_back - 0.5*dsigmaV_back);
    df_dsigmaV_back(i) = (f_plus-f_minus)/FD;
    % 2nd derivative
    if secondDerivative
        [df_dsigmaV_plus,~,~,~,~,df_dsigmaV_back_plus] = get_df_PlaneStressHardMixVAF_FD(theta,H_isotropic,sigmaV,gamma,sigmaV_back + 0.5*dsigmaV_back,FD,false);
        [df_dsigmaV_minus,~,~,~,~,df_dsigmaV_back_minus] = get_df_PlaneStressHardMixVAF_FD(theta,H_isotropic,sigmaV,gamma,sigmaV_back - 0.5*dsigmaV_back,FD,false);
        d2f_d2sigmaV_back(:,i) = (df_dsigmaV_back_plus - df_dsigmaV_back_minus)/FD;
        % mixed derivative
        d2f_dsigmaVdsigmaV_back(:,i) = (df_dsigmaV_plus - df_dsigmaV_minus)/FD;        
    end
    
end

end