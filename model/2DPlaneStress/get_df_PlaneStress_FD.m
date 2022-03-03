function [df_dsigmaV, d2f_d2sigmaV] = get_df_PlaneStress_FD(theta,sigmaV,FD,secondDerivative)
% get_df_PlaneStress_FD evaluates the yield function and its derivatives 
% using finite differences
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
% `FD` (_double_) - finite differences stepsize
%
% `secondDerivative` (_double_) - indicates whether second derivative
% should be computed
% 
% ## Output Arguments
% 
% `f` (_double_) - value of the yield function
% 
% `df_dsigmaV` (_double_) - first derivative of the yield function
% 
% `d2f_d2sigmaV` (_double_) - second derivative of the yield function

if nargin == 3
    secondDerivative = true;
end

df_dsigmaV = NaN(3,1);
d2f_d2sigmaV = NaN(3,3);

for i = 1:3
    dsigmaV = zeros(3,1);
    dsigmaV(i) = FD;
    
    % 1st derivative
    f_plus = get_f_PlaneStress(theta,sigmaV + 0.5*dsigmaV);
    f_minus = get_f_PlaneStress(theta,sigmaV - 0.5*dsigmaV);
    df_dsigmaV(i) = (f_plus-f_minus)/FD;
    % 2nd derivative
    if secondDerivative
        [df_dsigmaV_plus,~] = get_df_PlaneStress_FD(theta,sigmaV + 0.5*dsigmaV,FD,false);
        [df_dsigmaV_minus,~] = get_df_PlaneStress_FD(theta,sigmaV - 0.5*dsigmaV,FD,false);
        d2f_d2sigmaV(:,i) = (df_dsigmaV_plus - df_dsigmaV_minus)/FD;
    end
    
end

end