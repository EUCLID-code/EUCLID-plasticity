function [r, alpha] = get_r_alpha_PlaneStress(sigmaV)
% get_r_alpha_PlaneStress evaluates the Lode radius and Lode angle
%
% ## Comments
% 
% _none_
% 
% ## Input Arguments
% 
% `sigmaV` (_double_) - Cauchy stress under plane stress conditions in
% Voigt notation (sigma_11, sigma_22, sigma_22)
% 
% ## Output Arguments
% 
% `r` (_double_) - Lode radius
%
% `alpha` (_double_) - Lode angle
%

sigma = V2sigma(sigmaV);

%% Eigenvalues and Eigevectors
[~,eigenvalues] = eig_PlaneStress(sigma);

%% Lode Coordinates
R = [sqrt(2/3), -sqrt(1/6);...
    0, sqrt(1/2);];
eigenvalues_pi = R*eigenvalues;
r = sqrt(eigenvalues_pi(1)^2 + eigenvalues_pi(2)^2);
% check: norm_sigma_dev = norm([sigma_dev zeros(2,1);zeros(1,2) -1/3*sum(diag(sigma))],'fro'); 
alpha = atan2(eigenvalues_pi(2),eigenvalues_pi(1));
