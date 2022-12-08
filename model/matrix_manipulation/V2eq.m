function sigma_eq = V2eq(sigmaV)
% V2eq calculates the equivalent stress from the Cauchy stress vector (Voigt)
%

if length(sigmaV) == 6
    sigma_eq = sqrt(3/2)*sqrt(sum((sigmaV(1:3)-sum(sigmaV(1:3))/3).^2) + 2*sum((sigmaV(4:end)).^2));
elseif length(sigmaV) == 3
    % plane stress
    sigma_eq = sqrt(3/2)*sqrt(sum(([sigmaV(1:2);0]-sum(sigmaV(1:2))/3).^2) + 2*sigmaV(3).^2);
elseif length(sigmaV) == 4
    % plane strain (sigma_33 is not zero in general)
    sigma_eq = sqrt(3/2)*sqrt(sum(([sigmaV(1:2);sigmaV(4)]-sum([sigmaV(1:2);sigmaV(4)])/3).^2) + 2*sigmaV(3).^2);
else
    error('Input must have either 6 components in 3D or 3/4 components in 2D.')
end

end