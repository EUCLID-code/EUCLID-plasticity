function IV_eng = engineering2ordinary_strainV(dim)
% engineering2ordinary_strainV provides a matrix that transforms engineering strains to ordinary strains in Voigt notation
%
% ## Comments
%
% It is `epsilonV_ordinary = IV_eng * epsilonV_engineering`.
% `epsilonV_ordinary` is a vector containing the strain components in Voigt
% notation without a factor of two in front of the shear components.
% `epsilonV_engineering` is a vector containing the strain components in
% Voigt notation with a factor of two in front of the shear components.
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
% `IV_eng` (_double_) - matrix such that `epsilonV_ordinary = IV_eng * epsilonV_engineering`
% 

if nargin == 0
    dim = 3;
end

IV = identityV(dim);
IV_eng = diag((IV+1)/2);

end