function material = lib_VEVPHardMixLin_VEVPHardMixLin

% VE --> viscoelastic
% --> one Maxwell element
% VP --> viscoplastic
% --> linear viscoplastic evolution
% HardMixLin --> linear mixed hardening
% --> linear isotropic hardening
% --> linear kinematic hardening
material.name = 'lib_VEVPHardMixLin_VEVPHardMixLin';

%% Viscoelastic Parameters
material.G0 = 0.6; % initial shear modulus [GPa = kN/mm^2]
material.K0 = 1.3; % initial bulk modulus [GPa = kN/mm^2]
material.n_Maxwell_G = 1;
material.n_Maxwell_K = 1;
material.Gi = [0.35]; % [GPa = kN/mm^2]
material.gi = [110]; % [s]
material.Ki = [0.4]; % [GPa = kN/mm^2]
material.ki = [15]; % [s]
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
material.H_kin = 0.01; % [GPa = kN/mm^2]
% viscoplastic function
material.eta = 0.04; % [GPa s = kN/mm^2 s]
material.sigma_0 = 0.03; % yield stress [GPa = kN/mm^2]