function [f, df_dsigmaV, d2f_d2sigmaV] = get_df_PlaneStress_Tresca(sigma_0,sigmaV,ID)
% get_df_PlaneStress_Tresca evaluates the yield function and its derivatives
% for the Tresca model using finite differences
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

sigmaV = [sigmaV(1);sigmaV(2);0;0;0;sigmaV(3)];

%% Voigt
V = [ ...
    1, 1;
    2, 2;
    3, 3;
    2, 3;
    1, 3;
    1, 2;
    ];
I_V = [ ... % used to transform derivatives from the general notation to Voigt notation
    1, 0, 0, 0, 0, 0; ...
    0, 1, 0, 0, 0, 0; ...
    0, 0, 1, 0, 0, 0; ...
    0, 0, 0, 2, 0, 0; ...
    0, 0, 0, 0, 2, 0; ...
    0, 0, 0, 0, 0, 2; ...
    ];
sigma = V2sigma(sigmaV);

%% Eigenvalues and Eigevectors
[eigenvectors,eigenvalues] = eig(sigma);
% eigenvector1 = eigenvectors(:,1); eigenvector2 = eigenvectors(:,2); eigenvector3 = eigenvectors(:,3);
eigenvalues = diag(eigenvalues);
% change order of the eigenvalues and eigenvectors
sigma1 = eigenvalues(3); sigma2 = eigenvalues(2); sigma3 = eigenvalues(1);
eigenvalues(1) = sigma1; eigenvalues(2) = sigma2; eigenvalues(3) = sigma3;
e1 = eigenvectors(:,3); e2 = eigenvectors(:,2); e3 = eigenvectors(:,1);
eigenvectors(:,1) = e1; eigenvectors(:,2) = e2; eigenvectors(:,3) = e3;

%% Derivatives: Eigenvalues and Eigevectors
deigenvalues_dsigmaV = zeros(6,3); % 6x1 vector for every eigenvalue
for idxEV = 1:3 % loop over the three eigenvalues
    deigenvalues_dsigmaV(:,idxEV) = I_V*sigma2V(eigenvectors(:,idxEV)*eigenvectors(:,idxEV)');
end

% ATTENTION: symmetry w.r.t. idx1 and idx2 is enforced
deigenvectors_dsigma = zeros(3,3,3,3);
for idxEV = 1:3 % loop over the three eigenvectors
    for idxEVE = 1:3 % loop over the three entries of the eigenvectors
        for idx1 = 1:3 % loop over the entries of the stress tensor
            for idx2 = 1:3 % loop over the entries of the stress tensor
                for idxSum = 1:3 % sum
                    denominator = eigenvalues(idxEV) - eigenvalues(idxSum);
                    if abs(denominator) > 0
                        eigenvector_idxEV = eigenvectors(:,idxEV);
                        eigenvector_idxSum = eigenvectors(:,idxSum);
                        deigenvectors_dsigma(idxEV,idxEVE,idx1,idx2) = deigenvectors_dsigma(idxEV,idxEVE,idx1,idx2) + ...
                            0.5*eigenvector_idxSum(idxEVE)*eigenvector_idxSum(idx1)*eigenvector_idxEV(idx2)/denominator + ...
                            0.5*eigenvector_idxSum(idxEVE)*eigenvector_idxSum(idx2)*eigenvector_idxEV(idx1)/denominator;
                    end
                end
            end
        end
    end
end

% ATTENTION: derivatives must be transformed from the general notation to Voigt notation
d2eigenvalues_d2sigmaV = zeros(6,6,3); % 6x6 matrix for every eigenvalue
for idxEV = 1:3 % loop over the three eigenvalues
    for idxV1 = 1:6 % loop over the entries of the Voigt stress tensor
        for idxV2 = 1:6 % loop over the entries of the Voigt stress tensor
            V1 = I_V(idxV1,idxV1);
            V2 = I_V(idxV2,idxV2);
            idx1 = V(idxV1,1); idx2 = V(idxV1,2);
            idx3 = V(idxV2,1); idx4 = V(idxV2,2);
            eigenvector_idxEV = eigenvectors(:,idxEV);
            %     d2eigenvalues_d2sigmaV(idxV1,idxV2,idxEV) = V1*V2*d2eigenvalues_d2sigma(idxEV,idxS1,idxS2,idxS3,idxS4);
            d2eigenvalues_d2sigmaV(idxV1,idxV2,idxEV) = V1*V2* (...
                deigenvectors_dsigma(idxEV,idx1,idx3,idx4)*eigenvector_idxEV(idx2) + ...
                deigenvectors_dsigma(idxEV,idx2,idx3,idx4)*eigenvector_idxEV(idx1) );
        end
    end
end


%% Yield Function
if ID == 1
    f = sigma1 - sigma3 - sigma_0;
    df_dsigmaV = deigenvalues_dsigmaV(:,1) - deigenvalues_dsigmaV(:,3);
    d2f_d2sigmaV = squeeze(d2eigenvalues_d2sigmaV(:,:,1) - d2eigenvalues_d2sigmaV(:,:,3));
elseif ID == 2
    f = sigma2 - sigma3 - sigma_0;
    df_dsigmaV = deigenvalues_dsigmaV(:,2) - deigenvalues_dsigmaV(:,3);
    d2f_d2sigmaV = squeeze(d2eigenvalues_d2sigmaV(:,:,2) - d2eigenvalues_d2sigmaV(:,:,3));
elseif ID == 6
    f = sigma1 - sigma2 - sigma_0;
    df_dsigmaV = deigenvalues_dsigmaV(:,1) - deigenvalues_dsigmaV(:,2);
    d2f_d2sigmaV = squeeze(d2eigenvalues_d2sigmaV(:,:,1) - d2eigenvalues_d2sigmaV(:,:,2));
end

%% 2D
in_plane = [1,2,6];
df_dsigmaV = df_dsigmaV(in_plane);
d2f_d2sigmaV = d2f_d2sigmaV(in_plane,in_plane);


end
