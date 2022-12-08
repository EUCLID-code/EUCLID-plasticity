function [sigmaV,CC,alphaV,epsilonVvp,gamma,viscoelastic,converged_local] = local_problem_PlaneStrain_VEVPHardMixLin( ...
    Gi, ...
    gi, ...
    Ki, ...
    ki, ...
    Ginf, ...
    Kinf, ...
    H_iso, ... % isotropic hardening parameter
    H_kin, ... % kinematic hardening parameter
    eta, ...
    sigma_0, ...
    time_inc, ...
    epsilonV_2DPlaneStrain, ...
    alphaV_prev, ... % viscoelastic strain @ n
    epsilonVvp_prev, ... % viscoplastic strain @ n
    gamma_prev ... % plastic multiplier @ n
    )
% local_problem_PlaneStrain_VEVPHardMixLin solves the local problem
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
% load step in Voigt notation (epsilon_11, epsilon_22, 2*epsilon_12)
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
% notation (sigma_11, sigma_22, sigma_12)
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

%% Plane Strain Constraints
in_plane_components = [1,2,6];
epsilonV = zeros(6,1);
epsilonV(in_plane_components) = epsilonV_2DPlaneStrain;

[sigmaV,CC,alphaV,epsilonVvp,gamma,viscoelastic,converged_local] = local_problem_VEVPHardMixLin( ...
    Gi, ...
    gi, ...
    Ki, ...
    ki, ...
    Ginf, ...
    Kinf, ...
    H_iso, ... % isotropic hardening parameter
    H_kin, ... % kinematic hardening parameter
    eta, ...
    sigma_0, ...
    time_inc, ...
    epsilonV, ... % strain increment
    alphaV_prev, ... % viscoelastic strain @ n
    epsilonVvp_prev, ... % viscoplastic strain @ n
    gamma_prev ... % plastic multiplier @ n
    );

sigmaV = sigmaV(in_plane_components);
CC = CC(in_plane_components,in_plane_components);

end






















