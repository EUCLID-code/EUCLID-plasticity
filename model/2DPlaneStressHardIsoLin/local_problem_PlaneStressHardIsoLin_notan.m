function [sigmaV,epsilonVp,gamma,converged_local] = local_problem_PlaneStressHardIsoLin_notan( ...
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

converged_local = false;

%% Elastic Predictor
epsilonVp = epsilonVp_prev;
gamma = gamma_prev;
sigmaV = CPlaneStress * (epsilonV(1:3) - epsilonVp(1:3));
f = get_f_PlaneStressHardIsoLin(theta,H_isotropic,sigmaV,gamma);

%% Plastic Corrector
if f <= 0
    converged_local = true;
else
    Dgamma = 0;
    for k = 1:n_NR_local
        [f, df_dsigmaV, d2f_d2sigmaV, df_dgamma, d2f_d2gamma, d2f_dsigmaVdgamma] = get_df_PlaneStressHardIsoLin(theta,H_isotropic,sigmaV,gamma);
        R = - epsilonVp(1:3) + epsilonVp_prev(1:3) + Dgamma*df_dsigmaV;
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
%             jacobian = [ ...
%                 df_dsigmaV', df_dgamma; ...
%                 SPlaneStress+Dgamma*d2f_d2sigmaV, df_dsigmaV + Dgamma*d2f_dsigmaVdgamma];
%             rhs = [ ...
%                 f; ...
%                 R];
%             DfR = -jacobian\rhs;
%             sigmaV = sigmaV + DfR(1:end-1);
%             Dgamma = Dgamma + DfR(end);
%             gamma = gamma_prev + Dgamma;
%             epsilonVp(1:3) = - SPlaneStress*sigmaV + epsilonV(1:3);
            
            % Newton-Raphson (Simo & Hughes' version)
            A = inv(SPlaneStress + Dgamma*d2f_d2sigmaV); % saving costs by computing inverse once
            DDgamma = (f - df_dsigmaV'*A*R) / (df_dsigmaV'*A*df_dsigmaV);
            DepsilonVp = SPlaneStress*A*(R + DDgamma*df_dsigmaV);
            epsilonVp(1:3) = epsilonVp(1:3) + DepsilonVp;
            Dgamma = Dgamma + DDgamma;
            gamma = gamma_prev + Dgamma;
            sigmaV = CPlaneStress * (epsilonV(1:3) - epsilonVp(1:3));
            
        end
    end
    if ~converged_local
        disp('Newton-Raphson of the local problem did not converge.')
    end
end

end