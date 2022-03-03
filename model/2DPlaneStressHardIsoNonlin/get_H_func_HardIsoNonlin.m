function [H_func, dHdgamma, d2Hd2gamma] = get_H_func_HardIsoNonlin(H_isotropic,gamma)
% get_H_func_HardIsoNonlin computes the nonlinear hardening function ans
% its derivatives
%
% ## Comments
% 
% _none_
% 
% ## Input Arguments
% 
% `H_isotropic` (_double_) - hardening parameters
% 
% `gamma` (_double_) - plastic multiplier
% 
% ## Output Arguments
% 
% `H_func` (_double_) - hardening function
%
% `dHdgamma` (_double_) - hardening function derivative
%
% `d2Hd2gamma` (_double_) - hardening function 2nd derivative
%

H_func = H_isotropic(1)*gamma + H_isotropic(2)*(1-exp(-H_isotropic(3)*gamma));
dHdgamma = H_isotropic(1) + H_isotropic(2)*H_isotropic(3)*exp(-H_isotropic(3)*gamma);
d2Hd2gamma = -H_isotropic(2)*H_isotropic(3)^2*exp(-H_isotropic(3)*gamma);

% linear isotropic hardening
% H_func = H_isotropic(1)*gamma;
% dHdgamma = H_isotropic(1);
% d2Hd2gamma = 0;

end