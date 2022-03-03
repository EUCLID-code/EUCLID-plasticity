function [df_dsigmaV, d2f_d2sigmaV, df_dsigmaV_back, d2f_d2sigmaV_back, d2f_dsigmaVdsigmaV_back] = get_df_PlaneStressHardKinLin_FD(theta,sigmaV,sigmaV_back,FD,secondDerivative)

if nargin == 4
    secondDerivative = true;
end

df_dsigmaV = NaN(3,1);
d2f_d2sigmaV = NaN(3,3);
df_dsigmaV_back = NaN(3,1);
d2f_d2sigmaV_back = NaN(3,3);
d2f_dsigmaVdsigmaV_back = NaN(3,3);

% derivatives w.r.t. sigmaV
for i = 1:3
    dsigmaV = zeros(3,1);
    dsigmaV(i) = FD;
    
    % 1st derivative
    f_plus = get_f_PlaneStressHardKinLin(theta,sigmaV + 0.5*dsigmaV,sigmaV_back);
    f_minus = get_f_PlaneStressHardKinLin(theta,sigmaV - 0.5*dsigmaV,sigmaV_back);
    df_dsigmaV(i) = (f_plus-f_minus)/FD;
    % 2nd derivative
    if secondDerivative
        [df_dsigmaV_plus,~,~,~] = get_df_PlaneStressHardKinLin_FD(theta,sigmaV + 0.5*dsigmaV,sigmaV_back,FD,false);
        [df_dsigmaV_minus,~,~,~] = get_df_PlaneStressHardKinLin_FD(theta,sigmaV - 0.5*dsigmaV,sigmaV_back,FD,false);
        d2f_d2sigmaV(:,i) = (df_dsigmaV_plus - df_dsigmaV_minus)/FD;
    end
    
end

% derivatives w.r.t. sigmaV_back and mixed derivative
for i = 1:3
    dsigmaV_back = zeros(3,1);
    dsigmaV_back(i) = FD;
    
    % 1st derivative
    f_plus = get_f_PlaneStressHardKinLin(theta,sigmaV,sigmaV_back + 0.5*dsigmaV_back);
    f_minus = get_f_PlaneStressHardKinLin(theta,sigmaV,sigmaV_back - 0.5*dsigmaV_back);
    df_dsigmaV_back(i) = (f_plus-f_minus)/FD;
    % 2nd derivative
    if secondDerivative
        [df_dsigmaV_plus,~,df_dsigmaV_back_plus,~] = get_df_PlaneStressHardKinLin_FD(theta,sigmaV,sigmaV_back + 0.5*dsigmaV_back,FD,false);
        [df_dsigmaV_minus,~,df_dsigmaV_back_minus,~] = get_df_PlaneStressHardKinLin_FD(theta,sigmaV,sigmaV_back - 0.5*dsigmaV_back,FD,false);
        d2f_d2sigmaV_back(:,i) = (df_dsigmaV_back_plus - df_dsigmaV_back_minus)/FD;
        % mixed derivative
        d2f_dsigmaVdsigmaV_back(:,i) = (df_dsigmaV_plus - df_dsigmaV_minus)/FD;        
    end
    
end

end