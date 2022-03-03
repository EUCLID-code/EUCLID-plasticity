function [f, df_dsigmaV, d2f_d2sigmaV, df_dsigmaV_back, d2f_d2sigmaV_back, d2f_dsigmaVdsigmaV_back] = get_df_PlaneStressHardKinLin(theta,sigmaV,sigmaV_back)

sigmaV_rel = sigmaV - sigmaV_back; % relative stress = stress - back stress

[f, df_dsigmaV_rel, d2f_d2sigmaV_rel] = get_df_PlaneStress(theta,sigmaV_rel);

df_dsigmaV = df_dsigmaV_rel;
d2f_d2sigmaV = d2f_d2sigmaV_rel;
df_dsigmaV_back = - df_dsigmaV_rel;
d2f_d2sigmaV_back = d2f_d2sigmaV_rel;
d2f_dsigmaVdsigmaV_back = - d2f_d2sigmaV_rel;

end