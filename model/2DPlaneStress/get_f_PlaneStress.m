function f = get_f_PlaneStress(theta,sigmaV)
% get_f_PlaneStress evaluates the yield function
%
% ## Comments
% 
% _none_
% 
% ## Input Arguments
% 
% `theta` (_double_) - material parameters
% 
% `sigmaV` (_double_) - Cauchy stress under plane stress conditions in
% Voigt notation (sigma_11, sigma_22, sigma_22)
% 
% ## Output Arguments
% 
% `f` (_double_) - value of the yield function

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

%% Equivalent Stress
sigma_eq = sqrt(3/2)*r;

%% Yield Stress
ii = 3*(0:(length(theta)-1));
cc = cos(ii*alpha);
sigma_0 = sum(theta.*cc);

% sigma_0 = 0;
% for idx = 0:(length(theta)-1)
%     sigma_0 = sigma_0 + theta(idx+1)*cos(3*idx*alpha);
% end

%% Yield Function
f = sigma_eq - sigma_0;

end