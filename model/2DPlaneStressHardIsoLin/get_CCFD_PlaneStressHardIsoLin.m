function CC = get_CCFD_PlaneStressHardIsoLin( ...
    n_NR_local, ...
    tol_NR_local, ...
    CPlaneStress, ...
    SPlaneStress, ...
    theta, ...
    H_isotropic, ...
    epsilonV, ...
    epsilonVp_prev, ...
    gamma_prev, ...
    FD ...
    )

CC = zeros(3,3);
for i = 1:3
    depsilonV = zeros(size(epsilonV));
    depsilonV(i) = FD;
    
    epsilonV_plus = epsilonV + 0.5*depsilonV;
    sigmaV_plus = local_problem_PlaneStressHardIsoLin( ...
        n_NR_local, ...
        tol_NR_local, ...
        CPlaneStress, ...
        SPlaneStress, ...
        theta, ...
        H_isotropic, ...
        epsilonV_plus, ...
        epsilonVp_prev, ...
        gamma_prev ...
        );
    epsilonV_minus = epsilonV - 0.5*depsilonV;
    sigmaV_minus = local_problem_PlaneStressHardIsoLin( ...
        n_NR_local, ...
        tol_NR_local, ...
        CPlaneStress, ...
        SPlaneStress, ...
        theta, ...
        H_isotropic, ...
        epsilonV_minus, ...
        epsilonVp_prev, ...
        gamma_prev ...
        );
    
    CC(:,i) = (sigmaV_plus - sigmaV_minus)/FD;
end

end