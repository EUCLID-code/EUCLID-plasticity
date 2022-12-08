function IV_ord = ordinary2engineering_strainV(dim)
% ordinary2engineering_strainV provides a matrix that transforms ordinary strains to engineering strains in Voigt notation
%
% ## Comments
%
% It is `epsilonV_engineering = IV_ord * epsilonV_ordinary`.
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
% `IV_ord` (_double_) - matrix such that `epsilonV_engineering = IV_ord * epsilonV_ordinary`
% 

if nargin == 0
    dim = 3;
end

IV = identityV(dim);
IV_ord = diag((-IV+2));

end