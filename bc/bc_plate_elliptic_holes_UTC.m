function bc = bc_plate_elliptic_holes_UTC(mesh)

% ==================================================
% Uniaxial tension and compression.
% ==================================================

bc.name = 'plate_elliptic_holes_UTC';

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
    0.1; ... % case 1
    0.2; ... % case 2
    ];
bc.n_step = [ ...
    500; ... % case 1
    1000; ... % case 2
    ]; % number of load steps for each case

bc.displacement_increment = zeros(length(bc.max_displacement), length(bc.dof_displacement));
bc.displacement_increment(1,:) = [zeros(size(nodes_displacement)), ones(size(nodes_displacement))]; % case 1: vertical displacement, tension
bc.displacement_increment(2,:) = [zeros(size(nodes_displacement)), -ones(size(nodes_displacement))]; % case 2: vertical displacement, compression





