function bc = bc_rate_dependent_plate_elliptic_holes_UT_testcase(mesh)

% ==================================================
% Uniaxial tension.
% ==================================================

bc.name = 'bc_rate_dependent_plate_elliptic_holes_UT_testcase';

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
% rate-dependency requires to define the time interval during which the
% displacements or forces are applied
bc.max_time = [ ...
    1; ... % case 1
    1; ... % case 2
    4; ... % case 3
    ];
% the displacement max_displacement is applied during the time max_time
bc.max_displacement = [ ...
    6.25; ... % case 1
    6.25; ... % case 2
    0; ... % case 3
    ];
% time discretization
bc.n_step = [ ...
    5; ... % case 1
    5; ... % case 2
    20; ... % case 3
    ]; % number of load steps for each case

bc.displacement_increment = zeros(length(bc.max_displacement), length(bc.dof_displacement));
bc.displacement_increment(1,:) = [zeros(size(nodes_displacement)), ones(size(nodes_displacement))]; % case 1: vertical displacement, tension
bc.displacement_increment(2,:) = [zeros(size(nodes_displacement)), ones(size(nodes_displacement))]; % case 2: vertical displacement, tension
bc.displacement_increment(3,:) = [zeros(size(nodes_displacement)), zeros(size(nodes_displacement))]; % case 3: no displacement
 



