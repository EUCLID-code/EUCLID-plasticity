function [sigmaV,epsilonVp,converged_local] = local_problem_PlaneStress_notan( ...
    n_NR_local, ...
    tol_NR_local, ...
    CPlaneStress, ...
    SPlaneStress, ...
    theta, ...
    epsilonV, ...
    epsilonVp_prev ...
    )
% local_problem_PlaneStress_notan solves the local problem
%
% ## Comments
% 
% This function solves the local problem, i.e., the elastic
% predictor-plastic corrector return mapping algorithm. Note that it is not
% required to compute the elasto-plastic consistent tangent modulus in the
% inverse problem.
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
% `epsilonV` (_double_) - infinitesimal strain at the current load step in
% Voigt notation (epsilon_11, epsilon_22, 2*epsilon_12)
% 
% `epsilonVp_prev` (_double_) - plastic component of the infinitesimal
% strain at the previous load step in Voigt notation (epsilon_p_11,
% epsilon_p_22, 2*epsilon_p_12, epsilon_p_33)
% 
% ## Output Arguments
% 
% `sigmaV` (_double_) - Cauchy stress at the current load step in Voigt
% notation (sigma_11, sigma_22, sigma_22)
% 
% `epsilonVp` (_double_) - plastic component of the infinitesimal
% strain at the current load step in Voigt notation (epsilon_p_11,
% epsilon_p_22, 2*epsilon_p_12, epsilon_p_33)
% 
% `converged_local` (_logical_) - indicates whether the local problem
% converged

converged_local = false;

%% Elastic Predictor
epsilonVp = epsilonVp_prev;
sigmaV = CPlaneStress * (epsilonV(1:3) - epsilonVp(1:3));
f = get_f_PlaneStress(theta,sigmaV);

%% Plastic Corrector
if f <= 0
    converged_local = true;
else
    Dgamma = 0;
    for k = 1:n_NR_local
        [f, df_dsigmaV, d2f_d2sigmaV] = get_df_PlaneStress(theta,sigmaV);
        R = - epsilonVp(1:3) + epsilonVp_prev(1:3) + Dgamma*df_dsigmaV;
        if abs(f) < tol_NR_local && norm(R) < tol_NR_local
            if Dgamma < 0
                disp('Plastic multiplier is inadmissible.')
                converged_local = false;
                break
            end
            % Plastic Strain
            epsilonVp(4) = -(epsilonVp(1)+epsilonVp(2));
            % Continuum Consistent Tangent
            % (only applicable for very small "load steps"/"pseudo time steps")
%             CC_continuum = CPlaneStress ...
%                 - ((CPlaneStress*df_dsigmaV)*(CPlaneStress*df_dsigmaV)') ...
%                 / (df_dsigmaV'*CPlaneStress*df_dsigmaV);
            % Algorithmic Consistent Tangent (Simo & Hughes' version)
%             A = inv(SPlaneStress + Dgamma*d2f_d2sigmaV);
%             ss = df_dsigmaV'*A*df_dsigmaV;
%             if ss < 0
%                 disp('Square root does not exist.')
%                 converged_local = false;
%                 break
%             end
%             N = (A*df_dsigmaV) / sqrt(ss);
%             CC = A - N*N';
            converged_local = true;
            break
        else
            %             % Newton-Raphson (derived version)
            %             jacobian = [ ...
            %                 df_dsigmaV', 0; ...
            %                 SPlaneStress+Dgamma*d2f_d2sigmaV, df_dsigmaV];
            %             rhs = [ ...
            %                 f; ...
            %                 R];
            %             DfR = -jacobian\rhs;
            %             sigmaV = sigmaV + DfR(1:end-1);
            %             Dgamma = Dgamma + DfR(end);
            %             epsilonVp(1:3) = - SPlaneStress*sigmaV + epsilonV(1:3);
            
            % Newton-Raphson (Simo & Hughes' version)
            A = inv(SPlaneStress + Dgamma*d2f_d2sigmaV); % saving costs by computing inverse once
            DDgamma = (f - df_dsigmaV'*A*R) / (df_dsigmaV'*A*df_dsigmaV);
            DepsilonVp = SPlaneStress*A*(R + DDgamma*df_dsigmaV);
            epsilonVp(1:3) = epsilonVp(1:3) + DepsilonVp;
            Dgamma = Dgamma + DDgamma;
            sigmaV = CPlaneStress * (epsilonV(1:3) - epsilonVp(1:3));
        end
    end
    if ~converged_local
        disp('Newton-Raphson of the local problem did not converge.')
    end
end

end