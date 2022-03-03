function bc = bc_plate_elliptic_holes_shear(mesh)

% ==================================================
% Uniaxial tension and compression.
% ==================================================

bc.name = 'plate_elliptic_holes_shear';

%% Nodes
nodes_fix = mesh.bot_edge; % fixed nodes
nodes_displacement = mesh.top_edge; % displaced nodes

%% Degrees of Freedom
bc.dof_fix = [2*nodes_fix-1, 2*nodes_fix]; % fixed dof
bc.dof_displacement = [2*nodes_displacement-1, 2*nodes_displacement]; % displaced dof
bc.dof_reaction = [ ...
    2*nodes_displacement-1; ... % dof of 1st reaction force
    2*nodes_displacement; ... % dof of 2nd reaction force
    ]; % displaced dof

%% Dirichlet Boundary Conditions
bc.max_displacement = [ ...
    0.3; ... % case 1
    0.3; ... % case 2
    0.6; ... % case 3
    ];
bc.n_step = 500 / 0.3 * bc.max_displacement; % number of load steps for each case

bc.displacement_increment = zeros(length(bc.max_displacement), length(bc.dof_displacement));
bc.displacement_increment(1,:) = [ones(size(nodes_displacement)), zeros(size(nodes_displacement))]; % case 1: shear
bc.displacement_increment(2,:) = [zeros(size(nodes_displacement)), ones(size(nodes_displacement))]; % case 2: vertical displacement, tension
bc.displacement_increment(3,:) = [zeros(size(nodes_displacement)), -ones(size(nodes_displacement))]; % case 3: vertical displacement, compression





