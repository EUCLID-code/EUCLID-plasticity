function sigmaV_dev = V2dev(sigmaV)
% V2dev calculates the deviatoric stress from the Cauchy stress vector (Voigt)
%

if length(sigmaV) == 6
    sigmaV_dev = sigmaV;
    sigmaV_dev(1:3) = sigmaV_dev(1:3) - sum(sigmaV(1:3))/3;
elseif length(sigmaV) == 3
    % plane stress 
    % the out-of-plane deviatoric stress is not zero in
    % general
    sigmaV_dev = [sigmaV; 0];
    sigmaV_dev([1:2,4]) = sigmaV_dev([1:2,4]) - sum(sigmaV(1:2))/3;
elseif length(sigmaV) == 4
    sigmaV_dev = sigmaV;
    sigmaV_dev([1:2,4]) = sigmaV_dev([1:2,4]) - sum(sigmaV([1:2,4]))/3;
else
    error('Input must have either 6 components in 3D or 3/4 components in 2D.')
end

end