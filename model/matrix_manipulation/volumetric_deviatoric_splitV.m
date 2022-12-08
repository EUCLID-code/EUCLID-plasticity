function [IV_vol,IV_dev] = volumetric_deviatoric_splitV(dim)
% volumetric_deviatoric_splitV applies the volumetric deviatoric split to the stress tensor in Voigt notation
%
% ## Comments
%
% It is `sigmaV_vol = IV_vol * sigmaV` and `sigmaV_dev = IV_dev * sigmaV`.
%
% Assumptions:
% - Voigt notation
% 
% ## Input Arguments
% 
% `dim` (_double_) - dimension (default: 3)
% 
% ## Output Arguments
% 
% `IV_vol` (_double_) - matrix such that `sigmaV_vol = IV_vol * sigmaV`
% 
% `IV_dev` (_double_) - matrix such that `sigmaV_dev = IV_dev * sigmaV`
% 

if nargin == 0
    dim = 3;
    % if dim is not equal to three, this function should be used with care:
    % make sure that the plane stress / plane strain conditions are
    % fulfilled
end

IV = identityV(dim);
IV_vol = 1/3 * (IV*IV');
IV_dev = eye(length(IV)) - IV_vol;

end