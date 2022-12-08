function IV = identityV(dim)
% identity transforms the identity tensor from matrix to vector (Voigt) representation
% 
% ## Comments
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
% `IV` (_double_) - identity tensor as a vector in Voigt notation
% 

if nargin == 0
    dim = 3;
end

IV = sigma2V(eye(dim));

end