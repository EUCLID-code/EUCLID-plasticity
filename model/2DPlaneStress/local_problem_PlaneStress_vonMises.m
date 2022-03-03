function [sigmaV,CC,epsilonVp,elastic,converged_local] = local_problem_PlaneStress_vonMises( ...
    n_NR_local, ...
    tol_NR_local, ...
    E, ...
    nu, ...
    mu, ...
    CPlaneStress, ...
    SPlaneStress, ...
    sigma_0, ...
    epsilonV, ...
    epsilonVp_prev ...
    )
% local_problem_PlaneStress_vonMises solves the local problem for the von
% Mises model
%
% ## Comments
% 
% This function solves the local problem, i.e., the elastic
% predictor-plastic corrector return mapping algorithm.
% (see Neto et al. page 379)
%

elastic = true;
converged_local = false;

%%
P = 1/3*[ ...
    2, -1, 0; ...
    -1, 2, 0; ...
    0, 0, 6; ...
    ];

%% Elastic Predictor
epsilonVp = epsilonVp_prev;
sigmaV = CPlaneStress * (epsilonV(1:3) - epsilonVp(1:3));
a_1 = (sigmaV(1)+sigmaV(2))^2;
a_2 = (sigmaV(2)-sigmaV(1))^2;
a_3 = sigmaV(3)^2;
xi = 1/6*a_1 + 1/2*a_2 + 2*a_3;
Psi = 1/2*xi - 1/3*sigma_0^2;
CC = CPlaneStress;

%% Plastic Corrector
if Psi <= 0
    converged_local = true;
else
    elastic = false;
    Dgamma = 0;
    for k = 1:n_NR_local
        xi = ...
            + a_1/(6*(1+E*Dgamma/(3*(1-nu)))^2) ...
            + (1/2*a_2+2*a_3)/(1+2*mu*Dgamma)^2 ...
            ;
        Psi = 1/2*xi - 1/3*sigma_0^2;
        if abs(Psi) < tol_NR_local
            if Dgamma < 0
                disp('Plastic multiplier is inadmissible.')
                converged_local = false;
                break
            end
            % Stress
            EE = inv(SPlaneStress + Dgamma*P);
            A = EE * SPlaneStress;
            sigmaV = A*sigmaV;
            a_1 = (sigmaV(1)+sigmaV(2))^2;
            a_2 = (sigmaV(2)-sigmaV(1))^2;
            a_3 = sigmaV(3)^2;
            xi = 1/6*a_1 + 1/2*a_2 + 2*a_3;
            Psi = 1/2*xi - 1/3*sigma_0^2;
            if Psi > tol_NR_local
                disp('Local problem did not converge.')
                converged_local = false;
                break
            end
            % Plastic Strain
            epsilonVp(1:3) = -SPlaneStress*sigmaV + epsilonV(1:3);
            epsilonVp(4) = -(epsilonVp(1)+epsilonVp(2));
            epsilonVp_vol = 1/3*(epsilonVp(1)+epsilonVp(2)+epsilonVp(4));
            if abs(epsilonVp_vol) > tol_NR_local
                disp('Volumetric plastic strain occurs.')
            end
            % Algorithmic Consistent Tangent
            n = EE*P*sigmaV;
            alpha = 1/(sigmaV'*P*n);
            CC = EE - alpha*(n*n');
            converged_local = true;
            break
        else
            dxi = ...
                - a_1/(9*(1+E*Dgamma/(3*(1-nu)))^3) * E/(1-nu) ...
                - 2*mu*(a_2+4*a_3)/(1+2*mu*Dgamma)^3 ...
                ;
            df = 1/2*dxi;
            DDgamma = - Psi/df;
            Dgamma = Dgamma + DDgamma;
        end
    end
    if converged_local == false
        disp('Newton-Raphson of the local problem did not converge.')
    end
end

end