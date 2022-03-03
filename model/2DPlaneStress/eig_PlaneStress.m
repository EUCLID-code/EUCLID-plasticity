function [eigenvectors,eigenvalues] = eig_PlaneStress(sigma)
% eig_PlaneStress calculates eigenvalues and eigenvectors
%
% ## Comments
% 
% Eigenvalues and eigenvectors of a 2x2 symmetric matrix are calculated.
% 
% ## Input Arguments
% 
% `sigma` (_double_) - Cauchy stress under plane stress conditions given as
% a 2x2 symmetric matrix
% 
% ## Output Arguments
% 
% `eigenvectors` (_double_) - matrix containing the eigenvectors (each
% column corresponds to one eigenvector)
% 
% `eigenvalues` (_double_) - vector containing the eigenvalues in
% increasing order

a = sigma(1,1);
b = sigma(1,2); % is equal to c
d = sigma(2,2);
T = a+d;
D = a*d-b^2;

%% Eigenvalues and Eigenvectors
if b == 0
    L1 = a;
    L2 = d;
    V1 = [1; 0];
    V2 = [0; 1];
else
    L1 = T/2 + sqrt(T^2/4-D);
    L2 = T/2 - sqrt(T^2/4-D);
    V1 = [L1-d; b]; V1 = V1/norm(V1);
    V2 = [L2-d; b]; V2 = V2/norm(V2);
    %     V1 = [b; L1-a]; V1 = V1/norm(V1);
    %     V2 = [b; L2-a]; V2 = V2/norm(V2);
end

%% Order
if L1 <= L2
    eigenvalues = [L1; L2];
    eigenvectors = [V1, V2];
else
    eigenvalues = [L2; L1];
    eigenvectors = [V2, V1];
end

