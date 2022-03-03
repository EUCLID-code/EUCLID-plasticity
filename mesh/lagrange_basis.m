function [N,dNdxi] = lagrange_basis(type,coord)
% lagrange_basis returns the shape functions and gradients with respect to the reference coordinates
%

if strcmp(type,'Q4')

    xi = coord(1); eta=coord(2);
    N = 1/4*[ (1-xi)*(1-eta);
        (1+xi)*(1-eta);
        (1+xi)*(1+eta);
        (1-xi)*(1+eta)];
    dNdxi = 1/4*[-(1-eta), -(1-xi);
        1-eta, -(1+xi);
        1+eta, 1+xi;
        -(1+eta), 1-xi];

else
    error(['Element ',type,' is not implemented'])
end

end

