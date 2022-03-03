function [sigmaV,epsilonVp,gamma,sigmaV_back,converged_local] = local_problem_PlaneStressHardKinLin_notan( ...
    n_NR_local, ...
    tol_NR_local, ...
    CPlaneStress, ...
    SPlaneStress, ...
    theta, ...
    H_isotropic, ...
    H_kinematic, ...
    epsilonV, ...
    epsilonVp_prev, ...
    gamma_prev, ...
    sigmaV_back_prev ...
    )

converged_local = false;

%% Elastic Predictor
epsilonVp = epsilonVp_prev;
gamma = gamma_prev;
sigmaV_back = sigmaV_back_prev;
sigmaV = CPlaneStress * (epsilonV(1:3) - epsilonVp(1:3));
f = get_f_PlaneStressHardKinLin(theta,sigmaV,sigmaV_back);

df_dgamma = 0;
d2f_d2gamma = zeros(3,1);
d2f_dsigmaVdgamma = zeros(3,1);

%% Plastic Corrector
if f <= 0
    converged_local = true;
else
    Dgamma = 0;
    DsigmaV_back = zeros(3,1);
    for k = 1:n_NR_local
        [f, df_dsigmaV, d2f_d2sigmaV, df_dsigmaV_back, d2f_d2sigmaV_back, d2f_dsigmaVdsigmaV_back] = get_df_PlaneStressHardKinLin(theta,sigmaV,sigmaV_back);
        R = [ ...
            - epsilonVp(1:3) + epsilonVp_prev(1:3) + Dgamma*df_dsigmaV; ...
            - DsigmaV_back + H_kinematic*Dgamma*df_dsigmaV; ...
            ];
        if abs(f) < tol_NR_local && norm(R) < tol_NR_local
            if Dgamma < 0
                disp('Plastic multiplier is inadmissible.')
                converged_local = false;
                break
            end
            % Plastic Strain
            epsilonVp(4) = -(epsilonVp(1)+epsilonVp(2));
            converged_local = true;
            break
        else
            % Newton-Raphson (derived version)
            jacobian = [ ...
                df_dsigmaV', df_dgamma, df_dsigmaV_back'; ...
                SPlaneStress+Dgamma*d2f_d2sigmaV, df_dsigmaV + Dgamma*d2f_dsigmaVdgamma, Dgamma*d2f_dsigmaVdsigmaV_back; ...
                H_kinematic*Dgamma*d2f_d2sigmaV, H_kinematic*df_dsigmaV + H_kinematic*Dgamma*d2f_dsigmaVdgamma, - eye(3) + H_kinematic*Dgamma*d2f_dsigmaVdsigmaV_back; ...
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