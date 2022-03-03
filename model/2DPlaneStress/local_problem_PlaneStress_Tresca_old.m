function [sigmaV,CC,epsilonVp,elastic,converged_local] = local_problem_PlaneStress_Tresca( ...
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

elastic = true;
converged_local = false;

%% In-Plane Components
in_plane = [1,2,6];

%% Elastic Predictor
epsilonVp_trial = epsilonVp_prev;
sigmaV_trial = CPlaneStress * (epsilonV(1:3) - epsilonVp_trial(1:3));
CC = CPlaneStress;

%% Eigenvalues and Eigevectors
[~,eigenvalues_2D_trial] = eig_PlaneStress(V2sigma(sigmaV_trial));
% change order of the eigenvalues and eigenvectors
sigma1_2D_trial = eigenvalues_2D_trial(2); sigma2_2D_trial = eigenvalues_2D_trial(1);

% identify zero eigenvalue position
if sigma2_2D_trial > 0
    position = 3;
elseif sigma1_2D_trial < 0
    position = 1;
else
    position = 2;
end

sigmaV_3D_trial = [sigmaV_trial(1);sigmaV_trial(2);0;0;0;sigmaV_trial(3)];
[eigenvalues_trial, ~] = get_eigenvalues(sigmaV_3D_trial);
eigenvectors_trial = get_eigenvectors(sigmaV_3D_trial);
% change order of the eigenvalues and eigenvectors
sigma1_trial = eigenvalues_trial(3); sigma2_trial = eigenvalues_trial(2); sigma3_trial = eigenvalues_trial(1);

%% Eigenvalues and Eigevectors of the Deviatoric Stress
sigma_trial = V2sigma(sigmaV_3D_trial);
sigma_vol_trial = 1/3*sum(diag(sigma_trial))*eye(3);
p_trial = sigma_vol_trial(1,1);
sigma_dev_trial = sigma_trial - sigma_vol_trial;
sigmaV_dev_trial = sigma2V(sigma_dev_trial);
[eigenvalues_dev_trial, ~] = get_eigenvalues(sigmaV_dev_trial);
eigenvectors_dev_trial = get_eigenvectors(sigmaV_dev_trial);
% change order of the eigenvalues and eigenvectors
sigma1_dev_trial = eigenvalues_dev_trial(3); sigma2_dev_trial = eigenvalues_dev_trial(2); sigma3_dev_trial = eigenvalues_dev_trial(1);
e1_dev_trial = eigenvectors_dev_trial(:,3); e2_dev_trial = eigenvectors_dev_trial(:,2); e3_dev_trial = eigenvectors_dev_trial(:,1);
ee1_trial = e1_dev_trial*e1_dev_trial'; ee2_trial = e2_dev_trial*e2_dev_trial'; ee3_trial = e3_dev_trial*e3_dev_trial';

N1 = ee1_trial - ee3_trial;

%% Tests
sigma_trial_test = (sigma1_dev_trial + p_trial)*ee1_trial ...
        + (sigma2_dev_trial + p_trial)*ee2_trial ...
        + (sigma3_dev_trial + p_trial)*ee3_trial;
sigmaV_3D_trial_test = sigma2V(sigma_trial_test);
sigmaV_trial_test = sigmaV_3D_trial_test(in_plane);
epsilonVe_trial_test = SPlaneStress*sigmaV_trial_test;
epsilonVp_trial_test = epsilonV - epsilonVe_trial_test;
epsilonVp_trial_test(4) = - (epsilonVp_trial_test(1) + epsilonVp_trial_test(2));

switch position
    case 1
        if sigma1_trial ~= 0; error(''); end
        if sigma1_dev_trial ~= -p_trial; error(''); end
    case 2
        if sigma2_trial ~= 0; error(''); end
        if sigma2_dev_trial ~= -p_trial; error(''); end
    case 3
        if sigma3_trial ~= 0; error(''); end
        if sigma3_dev_trial ~= -p_trial; error(''); end
end  

%% Yield Function
f_trial = sigma1_trial - sigma3_trial - sigma_0;

%% Plastic Corrector
if f_trial < tol_NR_local
    %     disp('Admissible Elastic Predictor')
    sigmaV = sigmaV_trial;
    epsilonVp = epsilonVp_trial;
    converged_local = true;
else
    elastic = false;
    %     disp('Plastic Corrector')
    %% Return to Main Plane
    f_trial_test1 = sigma1_trial - sigma3_trial - sigma_0;
    f_trial_test2 = sigma1_dev_trial - sigma3_dev_trial - sigma_0;
    switch position
        case 1 % CASE C
            disp('First eigenvalue is zero')
            f_trial_test3 = - p_trial - sigma3_dev_trial - sigma_0;
            f_trial_test4 = max([abs(sigma1_trial - sigma2_trial);abs(sigma2_trial - sigma3_trial);abs(sigma1_trial - sigma3_trial)]) - sigma_0;
            d = 2*mu;
            Dgamma = d\f_trial;
            sigma1_dev = sigma1_dev_trial; % enforce plane stress
            sigma2_dev = sigma2_dev_trial;
            sigma3_dev = sigma3_dev_trial + 2*mu*Dgamma;
        case 2 % CASE A
            disp('Second eigenvalue is zero')
            d = 4*mu;
            Dgamma = d\f_trial;
            sigma1_dev = sigma1_dev_trial - 2*mu*Dgamma;
            sigma2_dev = sigma2_dev_trial; % enforce plane stress
            sigma3_dev = sigma3_dev_trial + 2*mu*Dgamma;
        case 3 % CASE B
            disp('Third eigenvalue is zero')
            f_trial_test3 = sigma1_dev_trial + p_trial - sigma_0;
            f_trial_test4 = max([abs(sigma1_trial - sigma2_trial);abs(sigma2_trial - sigma3_trial);abs(sigma1_trial - sigma3_trial)]) - sigma_0;
            d = 2*mu;
            Dgamma = d\f_trial;
            sigma1_dev = sigma1_dev_trial - 2*mu*Dgamma;
            sigma2_dev = sigma2_dev_trial;
            sigma3_dev = sigma3_dev_trial; % enforce plane stress
    end
    
    %% Check Admissibility of Return to Main Plane
    if sigma1_dev >= sigma2_dev && sigma2_dev >= sigma3_dev
        disp('Return to the main plane is admissible')
    else
        disp('Return to the main plane is inadmissible')
        
        if sigma1_dev - 2*sigma2_dev + sigma3_dev > 0
            disp('Return to the right corner')
            f6_trial = sigma1_trial - sigma2_trial - sigma_0;
            RHS = [f_trial; f6_trial];
            switch position
                case 1 % CASE C
                    d = [2, 0; 0, 2]*mu;
                    Dgamma = d\RHS;
                    sigma1_dev = sigma1_dev_trial;
                    sigma2_dev = sigma2_dev_trial + 2*mu*Dgamma(2);
                    sigma3_dev = sigma3_dev_trial + 2*mu*Dgamma(1);
                case 2 % CASE A
                    d = [4, 2; 2, 2]*mu;
                    Dgamma = d\RHS;
                    sigma1_dev = sigma1_dev_trial - 2*mu*sum(Dgamma);
                    sigma2_dev = sigma2_dev_trial;
                    sigma3_dev = sigma3_dev_trial + 2*mu*Dgamma(1);
                case 3 % CASE B
                    d = [2, 2; 2, 4]*mu;
                    Dgamma = d\RHS;
                    sigma1_dev = sigma1_dev_trial - 2*mu*sum(Dgamma);
                    sigma2_dev = sigma2_dev_trial + 2*mu*Dgamma(2);
                    sigma3_dev = sigma3_dev_trial;
            end
        else
            disp('Return to the left corner')
            f2_trial = sigma2_trial - sigma3_trial - sigma_0;
            RHS = [f_trial; f2_trial];
            switch position
                case 1 % CASE C
                    d = [2, 2; 2, 4]*mu;
                    Dgamma = d\RHS;
                    sigma1_dev = sigma1_dev_trial;
                    sigma2_dev = sigma2_dev_trial - 2*mu*Dgamma(2);
                    sigma3_dev = sigma3_dev_trial + 2*mu*sum(Dgamma);
                case 2 % CASE A
                    d = [4, 2; 2, 2]*mu;
                    Dgamma = d\RHS;
                    sigma1_dev = sigma1_dev_trial - 2*mu*Dgamma(1);
                    sigma2_dev = sigma2_dev_trial;
                    sigma3_dev = sigma3_dev_trial + 2*mu*sum(Dgamma);
                case 3 % CASE B
                    d = [2, 0; 0, 2]*mu;
                    Dgamma = d\RHS;
                    sigma1_dev = sigma1_dev_trial - 2*mu*Dgamma(1);
                    sigma2_dev = sigma2_dev_trial - 2*mu*Dgamma(2);
                    sigma3_dev = sigma3_dev_trial;
            end
            
        end
        
    end
    
    %% Update
    sigma = (sigma1_dev + p_trial)*ee1_trial ...
        + (sigma2_dev + p_trial)*ee2_trial ...
        + (sigma3_dev + p_trial)*ee3_trial;
    sigmaV_3D = sigma2V(sigma);
    sigmaV = sigmaV_3D(in_plane);
    epsilonVe = SPlaneStress*sigmaV;
    epsilonVp = epsilonV - epsilonVe;
    epsilonVp(4) = - (epsilonVp(1) + epsilonVp(2));
            
end

end
























