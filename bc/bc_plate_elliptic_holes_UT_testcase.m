function bc = bc_plate_elliptic_holes_UT_testcase(mesh)

% ==================================================
% Uniaxial tension.
% ==================================================

bc.name = 'plate_elliptic_holes_UT_testcase';

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
    ];
bc.n_step = [ ...
    10; ... % case 1
    ]; % number of load steps for each case

bc.displacement_increment = zeros(length(bc.max_displacement), length(bc.dof_displacement));
bc.displacement_increment(1,:) = [zeros(size(nodes_displacement)), ones(size(nodes_displacement))]; % case 1: vertical displacement, tension





