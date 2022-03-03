function f = get_f_PlaneStressHardKinLin(theta,sigmaV,sigmaV_back)

sigmaV_rel = sigmaV - sigmaV_back; % relative stress = stress - back stress

f = get_f_PlaneStress(theta,sigmaV_rel);

end