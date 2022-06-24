function [sigmaV,CC,epsilonVp,elastic,gamma,sigmaV_back,converged_local] = local_problem_PlaneStressHardMixVAF(n_NR_local,tol_NR_local,CPlaneStress,SPlaneStress,theta,H_isotropic,H_kinematic,epsilonV,epsilonVp_prev,gamma_prev,sigmaV_back_prev)
% local_problem_PlaneStressHardMixVAF solves the local problem
%
% ## Comments
% 
% This function solves the local problem, i.e., the elastic
% predictor-plastic corrector return mapping algorithm. Assumptions: plane
% stress, isotropic hardening (Voce), kinematic hardening
% (Armstrong-Frederick).
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
% `H_kinematic` (_double_) - kinematic hardening parameters
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
% `sigmaV_back_prev` (_double_) - back stress at the previous load step
%
% ## Output Arguments
% 
% `sigmaV` (_double_) - Cauchy stress at the current load step in Voigt
% notation (sigma_11, sigma_22, sigma_12)
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
% `sigma_back` (_double_) - back stress at the current load step
%
% `converged_local` (_logical_) - indicates whether the local problem
% converged

elastic = true;
converged_local = false;
H_kinematic_lin = H_kinematic(1); % linear hardening constant
H_kinematic_AF = H_kinematic(2); % Armstrong-Frederick hardening constant

%% Elastic Predictor
epsilonVp = epsilonVp_prev;
gamma = gamma_prev;
sigmaV_back = sigmaV_back_prev;
sigmaV = CPlaneStress * (epsilonV(1:3) - epsilonVp(1:3));
CC = CPlaneStress;
f = get_f_PlaneStressHardMixVAF(theta,H_isotropic,sigmaV,gamma,sigmaV_back);

%% Plastic Corrector
if f <= 0
    converged_local = true;
else
    elastic = false;
    Dgamma = 0;
    DsigmaV_back = zeros(3,1);
    for k = 1:n_NR_local
        [f, df_dsigmaV, d2f_d2sigmaV, df_dgamma, d2f_d2gamma, d2f_dsigmaVdgamma, df_dsigmaV_back, d2f_d2sigmaV_back, d2f_dsigmaVdsigmaV_back] = get_df_PlaneStressHardMixVAF(theta,H_isotropic,sigmaV,gamma,sigmaV_back);
        R = [ ...
            - epsilonVp(1:3) + epsilonVp_prev(1:3) + Dgamma*df_dsigmaV; ...
            - (1 + Dgamma*H_kinematic_AF)*sigmaV_back + sigmaV_back_prev + H_kinematic_lin*Dgamma*df_dsigmaV; ...
            ];
        if abs(f) < tol_NR_local && norm(R) < tol_NR_local
            if Dgamma < 0
                disp('Plastic multiplier is inadmissible.')
                converged_local = false;
                break
            end
            % Plastic Strain
            epsilonVp(4) = -(epsilonVp(1)+epsilonVp(2));
            % Algorithmic Consistent Tangent
            H1 = inv(eye(3) - inv(1 + Dgamma*H_kinematic_AF)*Dgamma*H_kinematic_lin*d2f_dsigmaVdsigmaV_back);
            H2 = inv(1 + Dgamma*H_kinematic_AF)*Dgamma*H_kinematic_lin*d2f_d2sigmaV;
            H3 = - (1 + Dgamma*H_kinematic_AF)^(-2)*H_kinematic_AF*(sigmaV_back_prev + Dgamma*H_kinematic_lin*df_dsigmaV) ...
                + inv(1 + Dgamma*H_kinematic_AF)*(H_kinematic_lin*df_dsigmaV + Dgamma*H_kinematic_lin*d2f_dsigmaVdgamma);
            H4 = Dgamma*d2f_d2sigmaV + Dgamma*d2f_dsigmaVdsigmaV_back*H1*H2;
            H5 = df_dsigmaV + Dgamma*d2f_dsigmaVdgamma + Dgamma*d2f_dsigmaVdsigmaV_back*H1*H3;
            %             H6 = inv(df_dgamma + df_dsigmaV_back'*H1*H3) * (- df_dsigmaV' - df_dsigmaV_back'*H1*H2); % H6 is not defined for certain choices of hardening parameters
            H7 = (df_dsigmaV' + df_dsigmaV_back'*H1*H2)*inv(SPlaneStress + H4);
            H8 = df_dgamma + df_dsigmaV_back'*H1*H3;
            H9 = inv(H7*H5-H8);
            H10 = inv(SPlaneStress + H4);
            % there are two ways to compute CC
            %             CC = inv(SPlaneStress + H4 + H5*H6); % H6 is not defined for certain choices of hardening parameters
            CC = H10 - H10*H5*H9*H7;
            converged_local = true;
            break
        else
            % Newton-Raphson (derived version)
            jacobian = [ ...
                df_dsigmaV', df_dgamma, df_dsigmaV_back'; ...
                SPlaneStress+Dgamma*d2f_d2sigmaV, df_dsigmaV + Dgamma*d2f_dsigmaVdgamma, Dgamma*d2f_dsigmaVdsigmaV_back; ...
                H_kinematic_lin*Dgamma*d2f_d2sigmaV, - H_kinematic_AF*sigmaV_back + H_kinematic_lin*df_dsigmaV + H_kinematic_lin*Dgamma*d2f_dsigmaVdgamma, - (1 + Dgamma*H_kinematic_AF)*eye(3) + H_kinematic_lin*Dgamma*d2f_dsigmaVdsigmaV_back; ...
                ];
            rhs = [ ...
                f; ...
                R];
            DfR = -jacobian\rhs;
            sigmaV = sigmaV + DfR(1:3);
            Dgamma = Dgamma + DfR(4);
            gamma = gamma_prev + Dgamma;
            DsigmaV_back = DsigmaV_back + DfR(end-2:end);
            sigmaV_back = sigmaV_back_prev + DsigmaV_back;
            epsilonVp(1:3) = - SPlaneStress*sigmaV + epsilonV(1:3);
            
        end
    end
    if ~converged_local
        disp('Newton-Raphson of the local problem did not converge.')
    end
end

end