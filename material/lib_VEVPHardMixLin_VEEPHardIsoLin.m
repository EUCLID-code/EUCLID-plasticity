function material = lib_VEVPHardMixLin_VEEPHardIsoLin

% VE --> viscoelastic
% --> one Maxwell element (only shear)
% EP --> elastoplastic
% HardIsoLin --> linear isotropic hardening
material.name = 'lib_VEVPHardMixLin_VEEPHardIsoLin';

min_parameter = 1e-6;
max_parameter = 1e9;

%% Viscoelastic Parameters
material.G0 = 0.6; % initial shear modulus [GPa = kN/mm^2]
material.K0 = 1.3; % initial bulk modulus [GPa = kN/mm^2]
material.n_Maxwell_G = 1;
material.n_Maxwell_K = 1;
material.Gi = [0.35]; % [GPa = kN/mm^2]
material.gi = [110]; % [s]
material.Ki = [0]; % [GPa = kN/mm^2]
material.ki = [max_parameter]; % [s]
material.Ginf = material.G0 - sum(material.Gi);
material.Kinf = material.K0 - sum(material.Ki);
material.eta_Gi = material.Gi.*material.gi; % viscosities
material.eta_Ki = material.Ki.*material.ki; % viscosities

%% Relaxation Functions
material.Gfun = @(t) material.Ginf + sum(material.Gi.*exp(-t./material.gi));
material.Kfun = @(t) material.Kinf + sum(material.Ki.*exp(-t./material.ki));

%% Viscoplastic Parameters
% hardening function
material.H_iso = 0.03; % [GPa = kN/mm^2]
material.H_kin = 0; % [GPa = kN/mm^2]
% viscoplastic function
material.eta = min_parameter; % [GPa s = kN/mm^2 s]
material.sigma_0 = 0.03; % yield stress [GPa = kN/mm^2]