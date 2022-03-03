function [f, df_dsigmaV, d2f_d2sigmaV, df_dgamma, d2f_d2gamma, d2f_dsigmaVdgamma] = get_df_PlaneStressHardIsoNonlin(theta,H_isotropic,sigmaV,gamma)
% get_df_PlaneStressHardIsoNonlin evaluates the yield function and its
% derivatives
%
% ## Comments
% 
% Assumptions: plane stress, isotropic hardening (Voce).
% 
% ## Input Arguments
% 
% `theta` (_double_) - material parameters
% 
% `H_isotropic` (_double_) - isotropic hardening parameters
%
% `sigmaV` (_double_) - Cauchy stress under plane stress conditions in
% Voigt notation (sigma_11, sigma_22, sigma_22)
% 
% `gamma` (_double_) - plastic multiplier
%
% ## Output Arguments
% 
% yield function value and derivatives of the yield function
% 

V = [ ...
    1, 1;
    2, 2;
    1, 2;
    ];
I_V = [ ... % used to transform derivatives from the general notation to Voigt notation
    1, 0, 0; ...
    0, 1, 0; ...
    0, 0, 2; ...
    ]; 
sigma = V2sigma(sigmaV);

%% Eigenvalues and Eigevectors
[eigenvectors,eigenvalues] = eig_PlaneStress(sigma);

%% Volumetric and Deviatoric Stress
sigma_vol = 1/3*sum(diag(sigma))*eye(2); % not considering the out-of-plane components
sigma_dev = sigma - sigma_vol; % not considering the out-of-plane components
sigmaV_dev = sigma2V(sigma_dev);
I_dev = [ ...
    2/3, -1/3, 0; ...
    -1/3, 2/3, 0; ...
    0, 0, 1; ...
    ];

%% Lode Coordinates
R = [sqrt(2/3), -sqrt(1/6);...
    0, sqrt(1/2);];
eigenvalues_pi = R*eigenvalues;
r = sqrt(eigenvalues_pi(1)^2 + eigenvalues_pi(2)^2);
% check: norm_sigma_dev = norm([sigma_dev zeros(2,1);zeros(1,2) -1/3*sum(diag(sigma))],'fro'); 
alpha = atan2(eigenvalues_pi(2),eigenvalues_pi(1));

%% Equivalent Stress
% ATTENTION: derivatives must be transformed from the general notation to Voigt notation
sigma_eq = sqrt(3/2)*r;
dsigma_eq_dsigmaV = sqrt(3/2)*I_V*sigmaV_dev./r;
d2sigma_eq_d2sigmaV = sqrt(3/2)*(I_V*I_dev./r - ((I_V*sigmaV_dev) * (I_V*sigmaV_dev)')./r^3);

%% Derivatives: Eigenvalues and Eigevectors
deigenvalues_dsigmaV = zeros(3,2); % 3x1 vector for every eigenvalue
for idxEV = 1:2 % loop over the two eigenvalues
    deigenvalues_dsigmaV(:,idxEV) = I_V*sigma2V(eigenvectors(:,idxEV)*eigenvectors(:,idxEV)');
end

% ATTENTION: symmetry w.r.t. idx1 and idx2 is enforced
deigenvectors_dsigma = zeros(2,2,2,2);
for idxEV = 1:2 % loop over the two eigenvectors
for idxEVE = 1:2 % loop over the two entries of the eigenvectors
for idx1 = 1:2 % loop over the entries of the stress tensor
for idx2 = 1:2 % loop over the entries of the stress tensor
for idxSum = 1:2 % sum
    denominator = eigenvalues(idxEV) - eigenvalues(idxSum);
    if abs(denominator) > 0
        eigenvector_idxEV = eigenvectors(:,idxEV);
        eigenvector_idxSum = eigenvectors(:,idxSum);
        deigenvectors_dsigma(idxEV,idxEVE,idx1,idx2) = deigenvectors_dsigma(idxEV,idxEVE,idx1,idx2) + ...
            0.5*eigenvector_idxSum(idxEVE)*eigenvector_idxSum(idx1)*eigenvector_idxEV(idx2)/denominator + ...
            0.5*eigenvector_idxSum(idxEVE)*eigenvector_idxSum(idx2)*eigenvector_idxEV(idx1)/denominator;
    end
end
end
end
end
end

% ATTENTION: derivatives must be transformed from the general notation to Voigt notation
d2eigenvalues_d2sigmaV = zeros(3,3,2); % 3x3 matrix for every eigenvalue
for idxEV = 1:2 % loop over the two eigenvalues
for idxV1 = 1:3 % loop over the entries of the Voigt stress tensor
for idxV2 = 1:3 % loop over the entries of the Voigt stress tensor
    V1 = I_V(idxV1,idxV1);
    V2 = I_V(idxV2,idxV2);
    idx1 = V(idxV1,1); idx2 = V(idxV1,2);
    idx3 = V(idxV2,1); idx4 = V(idxV2,2);
    eigenvector_idxEV = eigenvectors(:,idxEV);
%     d2eigenvalues_d2sigmaV(idxV1,idxV2,idxEV) = V1*V2*d2eigenvalues_d2sigma(idxEV,idxS1,idxS2,idxS3,idxS4);
    d2eigenvalues_d2sigmaV(idxV1,idxV2,idxEV) = V1*V2* (...
        deigenvectors_dsigma(idxEV,idx1,idx3,idx4)*eigenvector_idxEV(idx2) + ...
        deigenvectors_dsigma(idxEV,idx2,idx3,idx4)*eigenvector_idxEV(idx1) );
end
end
end

%% Derivatives: Lode Angle
denominator = ...
    (eigenvalues(1) - eigenvalues(2))^2 + eigenvalues(1)^2 + eigenvalues(2)^2;
dalpha_deigenvalues = sqrt(3)/denominator * [...
    -eigenvalues(2); ...
    eigenvalues(1); ...
    ];
d2alpha_d2eigenvalue1 = 2/denominator * ((-2)*eigenvalues(1)+eigenvalues(2)) * dalpha_deigenvalues(1);
d2alpha_d2eigenvalue2 = 2/denominator * (eigenvalues(1)+(-2)*eigenvalues(2)) * dalpha_deigenvalues(2);
d2alpha_d2eigenvalue12 = 2/denominator * (eigenvalues(1)+eigenvalues(2)) * sqrt(3)/denominator * (eigenvalues(2)-eigenvalues(1));
d2alpha_d2eigenvalues = [...
    d2alpha_d2eigenvalue1, d2alpha_d2eigenvalue12; ...
    d2alpha_d2eigenvalue12, d2alpha_d2eigenvalue2; ...
    ];

dalpha_dsigmaV = (dalpha_deigenvalues'*deigenvalues_dsigmaV')';

d2alpha_d2sigmaV = zeros(3,3);
for idxEV = 1:2 % loop over the two eigenvalues
    d2alpha_d2sigmaV = d2alpha_d2sigmaV + ...
        dalpha_deigenvalues(idxEV)*squeeze(d2eigenvalues_d2sigmaV(:,:,idxEV));
end

% Sum of d2alpha_d2eigenvalues times (deigenvalues_dsigmaV outer product deigenvalues_dsigmaV)
help_matrix = zeros(3,3);
for idxEV1 = 1:2 % loop over the two eigenvalues
for idxEV2 = 1:2 % loop over the two eigenvalues
    help_matrix = help_matrix + ...
        d2alpha_d2eigenvalues(idxEV1,idxEV2)*(squeeze(deigenvalues_dsigmaV(:,idxEV2))*squeeze(deigenvalues_dsigmaV(:,idxEV1))');
end
end

%% Yield Stress
ii = 3*(0:(length(theta)-1));
cc = cos(ii*alpha);
ss = sin(ii*alpha);
sigma_0 = sum(theta.*cc);
dsigma_0_dalpha = sum(theta.*ii.*(-ss));
d2sigma_0_d2alpha = sum(theta.*ii.^2.*(-cc));

% sigma_0 = 0;
% dsigma_0_dalpha = 0;
% d2sigma_0_d2alpha = 0;
% for idx = 0:(length(theta)-1)
%     sigma_0 = sigma_0 + theta(idx+1)*cos(3*idx*alpha);
%     dsigma_0_dalpha = dsigma_0_dalpha - theta(idx+1)*3*idx*sin(3*idx*alpha);
%     d2sigma_0_d2alpha = d2sigma_0_d2alpha - theta(idx+1)*3^2*idx^2*cos(3*idx*alpha);
% end

dsigma_0_dsigmaV = dsigma_0_dalpha*dalpha_dsigmaV;
d2sigma_0_d2sigmaV = ...
    d2sigma_0_d2alpha*(dalpha_dsigmaV*dalpha_dsigmaV') + ...
    dsigma_0_dalpha*help_matrix + ...
    dsigma_0_dalpha*d2alpha_d2sigmaV;

%% Hardening
[H_func, dHdgamma, d2Hd2gamma] = get_H_func_HardIsoNonlin(H_isotropic,gamma);

%% Yield Function
f = sigma_eq - (1 + H_func)*sigma_0;
df_dsigmaV = dsigma_eq_dsigmaV - (1 + H_func)*dsigma_0_dsigmaV;
d2f_d2sigmaV = d2sigma_eq_d2sigmaV - (1 + H_func)*d2sigma_0_d2sigmaV;
df_dgamma = - dHdgamma*sigma_0;
d2f_d2gamma = - d2Hd2gamma*sigma_0;
d2f_dsigmaVdgamma = - dHdgamma*dsigma_0_dsigmaV;

end