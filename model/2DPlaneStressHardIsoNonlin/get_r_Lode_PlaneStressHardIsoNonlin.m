function r = get_r_Lode_PlaneStressHardIsoNonlin(theta,H_isotropic,alpha,gamma)
% get_r_Lode_PlaneStressHardIsoNonlin evaluates the Lode radius at which
% f=0
% 

ii = 3*(0:(length(theta)-1));
cc = cos(ii*alpha);
sigma_0 = sum(theta.*cc); % yield stress

% sigma_0 = 0;
% for idx = 0:(length(theta)-1)
%     sigma_0 = sigma_0 + theta(idx+1)*cos(3*idx*alpha);
% end

%% Hardening
[H_func, ~, ~] = get_H_func_HardIsoNonlin(H_isotropic,gamma);

%% Yield Function
r = sqrt(2/3)*(1 + H_func)*sigma_0;

end