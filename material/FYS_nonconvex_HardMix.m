function material = FYS_nonconvex_HardMix

material.name = 'FYS_nonconvex_HardMix';

%% Scalar Parameters
material.E = 210; % Young's modulus (Steel)
material.nu = 0.3; % Poisson's ratio (Steel)
material.lambda = material.E*material.nu / ((1+material.nu)*(1-2*material.nu));
material.mu = material.E / (2*(1+material.nu)); % shear modulus
material.K = material.E / (3*(1-2*material.nu)); % bulk modulus
material.sigma_0 = 0.24; % yield stress (Neto et al.) [GPa = kN/mm^2]
material.H_isotropic = [60, 2, 500];
material.H_kinematic = [175, 700];

%% Elastic Properties (3D)
material.C = [ ...
    material.lambda+2*material.mu, material.lambda, material.lambda, 0, 0, 0;...
    material.lambda, material.lambda+2*material.mu, material.lambda, 0, 0, 0;...
    material.lambda, material.lambda, material.lambda+2*material.mu, 0, 0, 0;...
    0, 0, 0, material.mu, 0, 0;...
    0, 0, 0, 0, material.mu, 0;...
    0, 0, 0, 0, 0, material.mu;...
    ];

material.S = [ ...
    1/material.E, -material.nu/material.E, -material.nu/material.E, 0, 0, 0;...
    -material.nu/material.E, 1/material.E, -material.nu/material.E, 0, 0, 0;...
    -material.nu/material.E, -material.nu/material.E, 1/material.E, 0, 0, 0;...
    0, 0, 0, 1/material.mu, 0, 0;...
    0, 0, 0, 0, 1/material.mu, 0;...
    0, 0, 0, 0, 0, 1/material.mu;...
    ];

%% Elastic Properties (Plane Strain)
material.CPlaneStrain = [ ...
    material.lambda+2*material.mu, material.lambda, 0;...
    material.lambda, material.lambda+2*material.mu, 0;...
    0, 0, material.mu;...
    ];
% material.CPlaneStrain = E/((1+nu)*(1-2*nu)) * [ ...
%     1-nu, nu, 0;...
%     nu, 1-nu, 0;...
%     0, 0, 1/2-nu;...
%     ];

material.SPlaneStrain = inv(material.CPlaneStrain);

%% Elastic Properties (Plane Stress)
material.CPlaneStress = material.E/(1-material.nu^2) * [ ...
    1, material.nu, 0;...
    material.nu, 1, 0;...
    0, 0, 1/2*(1-material.nu);...
    ];

material.SPlaneStress = [ ...
    1/material.E, -material.nu/material.E, 0;...
    -material.nu/material.E, 1/material.E, 0;...
    0, 0, 1/material.mu;...
    ];
% material.SPlaneStress = inv(material.CPlaneStress);

material.oop = -material.nu/(1-material.nu); % factor to get out of plane elastic strain
% material.oop = -material.C(3,1)/material.C(3,3); % factor to get out of plane elastic strain

%% Plastic Properties
a(1+0*3)=0.17;
a(1+1*3)=0.07;
material.fourier_coeff = [a; zeros(size(a))];
material.theta = material.fourier_coeff(1,1:3:end);
