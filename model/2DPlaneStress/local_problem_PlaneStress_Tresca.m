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
% local_problem_PlaneStress_Tresca solves the local problem for the Tresca
% model
%
% ## Comments
% 
% This function solves the local problem, i.e., the elastic
% predictor-plastic corrector return mapping algorithm.
% (see Neto et al. page 379)
%

elastic = true;
converged_local = false;

%% In-Plane Components
in_plane = [1,2,6];

%% Elastic Predictor
epsilonVp_trial = epsilonVp_prev;
sigmaV_trial = CPlaneStress * (epsilonV(1:3) - epsilonVp_trial(1:3));
CC = CPlaneStress;
[f1_trial, df1_dsigmaV, d2f1_d2sigmaV] = get_df_PlaneStress_Tresca(sigma_0,sigmaV_trial,1);
[f2_trial, df2_dsigmaV, d2f2_d2sigmaV] = get_df_PlaneStress_Tresca(sigma_0,sigmaV_trial,2);
[f6_trial, df6_dsigmaV, d2f6_d2sigmaV] = get_df_PlaneStress_Tresca(sigma_0,sigmaV_trial,6);
if f2_trial < f6_trial
    right = 1; left = 0;
    fa_trial = f1_trial; dfa_dsigmaV = df1_dsigmaV; d2fa_d2sigmaV = d2f1_d2sigmaV;
    fb_trial = f6_trial; dfb_dsigmaV = df6_dsigmaV; d2fb_d2sigmaV = d2f6_d2sigmaV;
else
    right = 0; left = 1;
    fa_trial = f1_trial; dfa_dsigmaV = df1_dsigmaV; d2fa_d2sigmaV = d2f1_d2sigmaV;
    fb_trial = f2_trial; dfb_dsigmaV = df2_dsigmaV; d2fb_d2sigmaV = d2f2_d2sigmaV;
end

% %% Eigenvalues and Eigevectors
% [~,eigenvalues_2D_trial] = eig_PlaneStress(V2sigma(sigmaV_trial));
% % change order of the eigenvalues and eigenvectors
% sigma1_2D_trial = eigenvalues_2D_trial(2); sigma2_2D_trial = eigenvalues_2D_trial(1);
%
% % identify zero eigenvalue position
% if sigma2_2D_trial > 0
%     position = 3;
% elseif sigma1_2D_trial < 0
%     position = 1;
% else
%     position = 2;
% end
%
% sigmaV_3D_trial = [sigmaV_trial(1);sigmaV_trial(2);0;0;0;sigmaV_trial(3)];
% [eigenvalues_trial, ~] = get_eigenvalues(sigmaV_3D_trial);
% eigenvectors_trial = get_eigenvectors(sigmaV_3D_trial);
% % change order of the eigenvalues and eigenvectors
% sigma1_trial = eigenvalues_trial(3); sigma2_trial = eigenvalues_trial(2); sigma3_trial = eigenvalues_trial(1);
%
% %% Eigenvalues and Eigevectors of the Deviatoric Stress
% sigma_trial = V2sigma(sigmaV_3D_trial);
% sigma_vol_trial = 1/3*sum(diag(sigma_trial))*eye(3);
% p_trial = sigma_vol_trial(1,1);
% sigma_dev_trial = sigma_trial - sigma_vol_trial;
% sigmaV_dev_trial = sigma2V(sigma_dev_trial);
% [eigenvalues_dev_trial, ~] = get_eigenvalues(sigmaV_dev_trial);
% eigenvectors_dev_trial = get_eigenvectors(sigmaV_dev_trial);
% % change order of the eigenvalues and eigenvectors
% sigma1_dev_trial = eigenvalues_dev_trial(3); sigma2_dev_trial = eigenvalues_dev_trial(2); sigma3_dev_trial = eigenvalues_dev_trial(1);
% e1_dev_trial = eigenvectors_dev_trial(:,3); e2_dev_trial = eigenvectors_dev_trial(:,2); e3_dev_trial = eigenvectors_dev_trial(:,1);
% ee1_trial = e1_dev_trial*e1_dev_trial'; ee2_trial = e2_dev_trial*e2_dev_trial'; ee3_trial = e3_dev_trial*e3_dev_trial';
%
% %% Tests
% sigma_trial_test = (sigma1_dev_trial + p_trial)*ee1_trial ...
%     + (sigma2_dev_trial + p_trial)*ee2_trial ...
%     + (sigma3_dev_trial + p_trial)*ee3_trial;
% sigmaV_3D_trial_test = sigma2V(sigma_trial_test);
% sigmaV_trial_test = sigmaV_3D_trial_test(in_plane);
% epsilonVe_trial_test = SPlaneStress*sigmaV_trial_test;
% epsilonVp_trial_test = epsilonV - epsilonVe_trial_test;
% epsilonVp_trial_test(4) = - (epsilonVp_trial_test(1) + epsilonVp_trial_test(2));
%
% switch position
%     case 1
%         if sigma1_trial ~= 0; error(''); end
%         if sigma1_dev_trial ~= -p_trial; error(''); end
%     case 2
%         if sigma2_trial ~= 0; error(''); end
%         if sigma2_dev_trial ~= -p_trial; error(''); end
%     case 3
%         if sigma3_trial ~= 0; error(''); end
%         if sigma3_dev_trial ~= -p_trial; error(''); end
% end
%
% %% Yield Function
% f_trial = sigma1_trial - sigma3_trial - sigma_0;

%% Plastic Corrector
if f1_trial < tol_NR_local
    %     disp('Admissible Elastic Predictor')
    sigmaV = sigmaV_trial;
    epsilonVp = epsilonVp_trial;
    converged_local = true;
else
    elastic = false;
    %     disp('Plastic Corrector')
    %% Return to Main Plane
    Dgamma = 0;
    jacobian = [ ...
        df1_dsigmaV', 0; ...
        SPlaneStress+Dgamma*d2f1_d2sigmaV, df1_dsigmaV];
    R = - epsilonVp_trial(1:3) + epsilonVp_prev(1:3) + Dgamma*df1_dsigmaV;
    rhs = [ ...
        f1_trial; ...
        R];
    DfR = -jacobian\rhs;
    Dgamma = Dgamma + DfR(end);
    epsilonVp_new = epsilonVp_trial(1:3) + Dgamma*df1_dsigmaV;
    
    %% Update
    epsilonVe_new = epsilonV - epsilonVp_new;
    sigmaV_new = CPlaneStress * epsilonVe_new;
    epsilonVp_new(4) =  - (epsilonVp_new(1) + epsilonVp_new(2));
    sigmaV = sigmaV_new;
    epsilonVp = epsilonVp_new;
    
    %% Eigenvalues and Eigevectors of the Deviatoric Stress
    sigmaV_3D_new = [sigmaV_new(1);sigmaV_new(2);0;0;0;sigmaV_new(3)];
    sigma_new = V2sigma(sigmaV_3D_new);
    sigma_vol_new = 1/3*sum(diag(sigma_new))*eye(3);
    p_new = sigma_vol_new(1,1);
    sigma_dev_new = sigma_new - sigma_vol_new;
    sigmaV_dev_new = sigma2V(sigma_dev_new);
    [eigenvalues_dev_new, ~] = get_eigenvalues(sigmaV_dev_new);
    % change order of the eigenvalues and eigenvectors
    sigma1_dev_new = eigenvalues_dev_new(3); sigma2_dev_new = eigenvalues_dev_new(2); sigma3_dev_new = eigenvalues_dev_new(1);
    f_trial_new = max([abs(sigma1_dev_new - sigma2_dev_new);abs(sigma2_dev_new - sigma3_dev_new);abs(sigma1_dev_new - sigma3_dev_new)]) - sigma_0;
    
    %% Check Admissibility of Return to Main Plane
    %     if sigma1_dev_new >= sigma2_dev_new && sigma2_dev_new >= sigma3_dev_new
    if f_trial_new < tol_NR_local
%         disp('Return to the main plane is admissible')
        % Algorithmic Consistent Tangent (Simo & Hughes' version)
%         A = inv(SPlaneStress + Dgamma*d2f1_d2sigmaV);
%         ss = df1_dsigmaV'*A*df1_dsigmaV;
%         if ss < 0
%             disp('Square root does not exist.')
%             converged_local = false;
%             return
%         end
%         N = (A*df1_dsigmaV) / sqrt(ss);
%         CC = A - N*N';
        converged_local = true;
    else
%         disp('Return to the main plane is inadmissible')
        %% Return to Corner
        Dgamma_a = 0;
        Dgamma_b = 0;
        jacobian = [ ...
            dfa_dsigmaV', 0, 0; ...
            dfb_dsigmaV', 0, 0; ...
            SPlaneStress+Dgamma_a*d2fa_d2sigmaV+Dgamma_b*d2fb_d2sigmaV, dfa_dsigmaV, dfb_dsigmaV];
        R = - epsilonVp_trial(1:3) + epsilonVp_prev(1:3) + Dgamma_a*dfa_dsigmaV + Dgamma_b*dfb_dsigmaV;
        rhs = [ ...
            fa_trial; ...
            fb_trial; ...
            R];
        DfR = -jacobian\rhs;
        Dgamma_a = Dgamma_a + DfR(end-1);
        Dgamma_b = Dgamma_b + DfR(end);
        epsilonVp_new = epsilonVp_trial(1:3) + Dgamma_a*dfa_dsigmaV + Dgamma_b*dfb_dsigmaV;
        
        %% Update
        epsilonVe_new = epsilonV - epsilonVp_new;
        sigmaV_new = CPlaneStress * epsilonVe_new;
        epsilonVp_new(4) =  - (epsilonVp_new(1) + epsilonVp_new(2));
        sigmaV = sigmaV_new;
        epsilonVp = epsilonVp_new;
        
        %% Eigenvalues and Eigevectors of the Deviatoric Stress
        sigmaV_3D_new = [sigmaV_new(1);sigmaV_new(2);0;0;0;sigmaV_new(3)];
        sigma_new = V2sigma(sigmaV_3D_new);
        sigma_vol_new = 1/3*sum(diag(sigma_new))*eye(3);
        p_new = sigma_vol_new(1,1);
        sigma_dev_new = sigma_new - sigma_vol_new;
        sigmaV_dev_new = sigma2V(sigma_dev_new);
        [eigenvalues_dev_new, ~] = get_eigenvalues(sigmaV_dev_new);
        % change order of the eigenvalues and eigenvectors
        sigma1_dev_new = eigenvalues_dev_new(3); sigma2_dev_new = eigenvalues_dev_new(2); sigma3_dev_new = eigenvalues_dev_new(1);
        f_trial_new = max([abs(sigma1_dev_new - sigma2_dev_new);abs(sigma2_dev_new - sigma3_dev_new);abs(sigma1_dev_new - sigma3_dev_new)]) - sigma_0;
        
        %% Check Admissibility of Return to Corner
        if f_trial_new < tol_NR_local
%             if right
%                 disp('Return to the right corner is admissible')
%             end
%             if left
%                 disp('Return to the left corner is admissible')
%             end
            converged_local = true;
        end
        
    end
    
end
























