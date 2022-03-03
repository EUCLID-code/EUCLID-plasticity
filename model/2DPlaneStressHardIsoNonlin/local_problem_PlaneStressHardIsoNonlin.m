function [sigmaV,CC,epsilonVp,elastic,gamma,converged_local] = local_problem_PlaneStressHardIsoNonlin( ...
    n_NR_local, ...
    tol_NR_local, ...
    CPlaneStress, ...
    SPlaneStress, ...
    theta, ...
    H_isotropic, ...
    epsilonV, ...
    epsilonVp_prev, ...
    gamma_prev ...
    )
% local_problem_PlaneStressHardIsoNonlin solves the local problem
%
% ## Comments
% 
% This function solves the local problem, i.e., the elastic
% predictor-plastic corrector return mapping algorithm. Assumptions: plane
% stress, isotropic hardening (Voce).
% 
% ## Input Arguments
% 
% `n_NR_local` (_double_) - maximum number of Newton-Raphson iterations
% used for solving the local problem (return mapping algorithm)
% 
% `tol_NR_local` (_double_) - stopping tolerance for the Newton-Raphson
% iteration used for solving the local problem (return mapping algorithm)
% 
% `CPlaneStress` (_double_) - elastic properties (stiffness matrix)
% 
% `SPlaneStress` (_double_) - elastic properties (compliance matrix)
% 
% `theta` (_double_) - material parameters
%
% `H_isotropic` (_double_) - isotropic hardening parameters
% 
% `epsilonV` (_double_) - infinitesimal strain at the current load step in
% Voigt notation (epsilon_11, epsilon_22, 2*epsilon_12)
% 
% `epsilonVp_prev` (_double_) - plastic component of the infinitesimal
% strain at the previous load step in Voigt notation (epsilon_p_11,
% epsilon_p_22, 2*epsilon_p_12, epsilon_p_33)
% 
% `gamma_prev` (_double_) - plastic multiplier at the previous load step
%
% ## Output Arguments
% 
% `sigmaV` (_double_) - Cauchy stress at the current load step in Voigt
% notation (sigma_11, sigma_22, sigma_22)
%
% `CC` (_double_) - consistent tangent modulus
% 
% `epsilonVp` (_double_) - plastic component of the infinitesimal
% strain at the current load step in Voigt notation (epsilon_p_11,
% epsilon_p_22, 2*epsilon_p_12, epsilon_p_33)
%
% `elastic` (_logical_) - indicates whether the load step was purely
% elastic
% 
% `gamma` (_double_) - plastic multiplier at the current load step
%
% `converged_local` (_logical_) - indicates whether the local problem
% converged

elastic = true;
converged_local = false;

%% Elastic Predictor
epsilonVp = epsilonVp_prev;
gamma = gamma_prev;
sigmaV = CPlaneStress * (epsilonV(1:3) - epsilonVp(1:3));
CC = CPlaneStress;
f = get_f_PlaneStressHardIsoNonlin(theta,H_isotropic,sigmaV,gamma);

%% Plastic Corrector
if f <= 0
    converged_local = true;
else
    elastic = false;
    Dgamma = 0;
    for k = 1:n_NR_local
        [f, df_dsigmaV, d2f_d2sigmaV, df_dgamma, d2f_d2gamma, d2f_dsigmaVdgamma] = get_df_PlaneStressHardIsoNonlin(theta,H_isotropic,sigmaV,gamma);
        R = - epsilonVp(1:3) + epsilonVp_prev(1:3) + Dgamma*df_dsigmaV;
        if abs(f) < tol_NR_local && norm(R) < tol_NR_local
            if Dgamma < 0
                disp('Plastic multiplier is inadmissible.')
                converged_local = false;
                break
            end
            % Plastic Strain
            epsilonVp(4) = -(epsilonVp(1)+epsilonVp(2));
            % Algorithmic Consistent Tangent
            A = inv(SPlaneStress + Dgamma*d2f_d2sigmaV);
            ss = df_dsigmaV'*A*df_dsigmaV + Dgamma*df_dsigmaV'*A*d2f_dsigmaVdgamma - df_dgamma;
            if ss < 0
                disp('Square root does not exist.')
                converged_local = false;
                break
            end
            N1 = (A*df_dsigmaV) / sqrt(ss);
            N2 = (A*(df_dsigmaV + Dgamma*d2f_dsigmaVdgamma)) / sqrt(ss);
            CC = A - N2*N1';
            converged_local = true;
            break
        else
            % Newton-Raphson (derived version)
            jacobian = [ ...
                df_dsigmaV', df_dgamma; ...
                SPlaneStress+Dgamma*d2f_d2sigmaV, df_dsigmaV + Dgamma*d2f_dsigmaVdgamma];
            rhs = [ ...
                f; ...
                R];
            DfR = -jacobian\rhs;
            sigmaV = sigmaV + DfR(1:end-1);
            Dgamma = Dgamma + DfR(end);
            gamma = gamma_prev + Dgamma;
            epsilonVp(1:3) = - SPlaneStress*sigmaV + epsilonV(1:3);
        end
    end
    if ~converged_local
        disp('Newton-Raphson of the local problem did not converge.')
    end
end

end