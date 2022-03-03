function f = get_f_Lode_PlaneStressHardIsoNonlin(theta,H_isotropic,r,alpha,gamma)

sigma_eq = sqrt(3/2)*r; % eqivalent stress

%% Yield Stress
ii = 3*(0:(length(theta)-1));
cc = cos(ii*alpha);
sigma_0 = sum(theta.*cc);

% sigma_0 = 0;
% for idx = 0:(length(theta)-1)
%     sigma_0 = sigma_0 + theta(idx+1)*cos(3*idx*alpha);
% end

%% Hardening
[H_func, ~, ~] = get_H_func_HardIsoNonlin(H_isotropic,gamma);

%% Yield Function
f = sigma_eq - (1 + H_func)*sigma_0;

end