function [sigmaV,CC,alphaV,epsilonVvp,gamma,viscoelastic,converged_local] = local_problem_VEVPHardMixLin(Gi,gi,Ki,ki,Ginf,Kinf,H_iso,H_kin,eta,sigma_0,time_inc,epsilonV,alphaV_prev,epsilonVvp_prev,gamma_prev)
% local_problem_VEVPHardMixLin solves the local problem
%
% ## Comments
%
% This function solves the local problem, i.e., the
% viscoelastic-viscoplastic constitutive equations with mixed isotropic and
% kinematic hardening.
%
% Note:
%
% `epsilon_vp = alpha_1 = alpha_3`
%
% `gamma = sqrt(3/2) * alpha_2`
% 
% ## Input Arguments
% 
% `Gi` (_double_) - viscoelastic material parameters
% 
% `gi` (_double_) - viscoelastic material parameters
% 
% `Ki` (_double_) - viscoelastic material parameters
% 
% `ki` (_double_) - viscoelastic material parameters
% 
% `Ginf` (_double_) - viscoelastic material parameter
%
% `Kinf` (_double_) - viscoelastic material parameter
%
% `H_iso` (_double_) - isotropic hardening parameter
% 
% `H_kin` (_double_) - kinematic hardening parameter
% 
% `eta` (_double_) - viscoplastic material parameter
%
% `sigma_0` (_double_) - yield stress
% 
% `time_inc` (_double_) - time increment
% 
% `epsilonV_2DPlaneStrain` (_double_) - infinitesimal strain at the current
% load step in Voigt notation
% 
% `alphaV_prev` (_double_) - viscoelastic internal variables at the
% previous load step
% 
% `epsilonVvp_prev` (_double_) - viscoplastic internal variables at the
% previous load step
%
% `gamma_prev` (_double_) - plastic multiplier at the previous load step
%
% ## Output Arguments
% 
% `sigmaV` (_double_) - Cauchy stress at the current load step in Voigt
% notation
%
% `CC` (_double_) - consistent tangent modulus
% 
% `alphaV` (_double_) - viscoelastic internal variables at the
% current load step
%
% `epsilonVvp` (_double_) - viscoplastic internal variables at the
% current load step
% 
% `gamma` (_double_) - plastic multiplier at the current load step
%
% `viscoelastic` (_logical_) - indicates whether the load step was purely
% viscoelastic
%
% `converged_local` (_logical_) - indicates whether the local problem
% converged

% caution: this implementation assumes gi > 0 and ki > 0
TOL = 0;

%% Volumetric-Deviatoric Split
[IV_vol,IV_dev] = volumetric_deviatoric_splitV(3);
IV_eng = engineering2ordinary_strainV(3);
IV_ord = ordinary2engineering_strainV(3);
epsilonVve = epsilonV - epsilonVvp_prev;
epsilonV_vol = IV_vol*epsilonVve;
epsilonV_dev = IV_dev*epsilonVve;

%% Useful Parameters
Gigi = time_inc./gi;
Kiki = time_inc./ki;
G_epsilon = 2*Ginf;
for idx = 1:length(Gi)
    if Gi(idx) > TOL
        G_epsilon = G_epsilon + 2*Gi(idx) * (1 - Gigi(idx)/(1 + Gigi(idx)));
    end
end

%% Update
alphaV = zeros(size(alphaV_prev));
sigmaV = 2*Ginf*IV_eng*epsilonV_dev + 3*Kinf*IV_eng*epsilonV_vol;

for idx = 1:length(Gi)
    if Gi(idx) > TOL
        alphaV_dev_prev = IV_dev*alphaV_prev;
        alphaV_dev = (alphaV_dev_prev + Gigi(idx) * epsilonV_dev) / (1 + Gigi(idx));
        alphaV(:,idx) = alphaV(:,idx) + alphaV_dev;
        sigmaV = sigmaV + 2*Gi(idx)*IV_eng*(epsilonV_dev - alphaV_dev);
    end
end
for idx = 1:length(Ki)
    if Ki(idx) > TOL
        alphaV_vol_prev = IV_vol*alphaV_prev;
        alphaV_vol = (alphaV_vol_prev + Kiki(idx) * epsilonV_vol) / (1 + Kiki(idx));
        alphaV(:,idx) = alphaV(:,idx) + alphaV_vol;
        sigmaV = sigmaV + 3*Ki(idx)*IV_eng*(epsilonV_vol - alphaV_vol);
    end
end

%% Trial
sigma_eq = V2eq(sigmaV - H_kin*epsilonVvp_prev); % relative equivalent stress
f_trial = sigma_eq - sigma_0 - sqrt(2/3)*H_iso*gamma_prev;
sigma_dev = V2dev(sigmaV - H_kin*epsilonVvp_prev); % relative deviatoric stress
sigma_dev_norm = sqrt(2/3)*sigma_eq; % norm of the relative deviatoric stress
N = sigma_dev / sigma_dev_norm; % director of the relative deviatoric stress

%% Consistent Tangent
G_sum = 2*Ginf;
for idx = 1:length(Gi)
    if Gi(idx) > TOL
        G_sum = G_sum + 2*Gi(idx) * (1 - Gigi(idx)/(1 + Gigi(idx)));
    end
end
K_sum = 3*Kinf;
for idx = 1:length(Ki)
    if Ki(idx) > TOL
        K_sum = K_sum + 3*Ki(idx) * (1 - Kiki(idx)/(1 + Kiki(idx)));
    end
end
CC_trial = G_sum*IV_eng*IV_dev + K_sum*IV_eng*IV_vol;

%% Update
if f_trial <= 0
    viscoelastic = true;
    converged_local = true;

    epsilonVvp = epsilonVvp_prev;
    gamma = gamma_prev;
    CC = CC_trial;
else
    viscoelastic = false;
    converged_local = true;

    f_denominator = sqrt(2/3)*eta/time_inc + sqrt(2/3)*H_iso + sqrt(3/2)*(G_epsilon + H_kin);
    gamma_inc = f_trial / f_denominator;
    gamma = gamma_prev + gamma_inc;
    epsilonVvp_inc = gamma_inc * IV_ord * N;
    epsilonVvp = epsilonVvp_prev + epsilonVvp_inc;
    epsilonVve = epsilonV - epsilonVvp;

    % viscoelastic stress computation
    epsilonV_vol = IV_vol*epsilonVve;
    epsilonV_dev = IV_dev*epsilonVve;
    alphaV = zeros(size(alphaV_prev));
    sigmaV = 2*Ginf*IV_eng*epsilonV_dev + 3*Kinf*IV_eng*epsilonV_vol;

    for idx = 1:length(Gi)
        if Gi(idx) > TOL
            alphaV_dev_prev = IV_dev*alphaV_prev;
            alphaV_dev = (alphaV_dev_prev + Gigi(idx) * epsilonV_dev) / (1 + Gigi(idx));
            alphaV(:,idx) = alphaV(:,idx) + alphaV_dev;
            sigmaV = sigmaV + 2*Gi(idx)*IV_eng*(epsilonV_dev - alphaV_dev);
        end
    end
    for idx = 1:length(Ki)
        if Ki(idx) > TOL
            alphaV_vol_prev = IV_vol*alphaV_prev;
            alphaV_vol = (alphaV_vol_prev + Kiki(idx) * epsilonV_vol) / (1 + Kiki(idx));
            alphaV(:,idx) = alphaV(:,idx) + alphaV_vol;
            sigmaV = sigmaV + 3*Ki(idx)*IV_eng*(epsilonV_vol - alphaV_vol);
        end
    end

    CC = K_sum*IV_eng*IV_vol ...
        + (1-G_epsilon*gamma_inc/sigma_dev_norm)*G_sum*IV_eng*IV_dev ...
        - (G_epsilon/f_denominator*(sigma_0 + sqrt(2/3)*H_iso*gamma_prev)/sigma_dev_norm^3)*IV_ord*CC_trial*(sigma_dev*sigma_dev');

end

end






















