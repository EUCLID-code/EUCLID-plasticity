function bc = bc_rate_dependent_plate_elliptic_holes_UTC_Nphases(mesh)

% ==================================================
% Uniaxial tension and compression.
% ==================================================

bc.name = 'bc_rate_dependent_plate_elliptic_holes_UTC_Nphases';

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
% displacements or forces are applied in n phases
% every phase consists of a loading period (during which the displacements
% or forces are applied) and a relaxation period (during which no
% displacements or forces are applied)
bc.n_phase = 20;
% rate-dependency requires to define the time interval during which the
% displacements or forces are applied
bc.max_time_load = logspace(-2,2,bc.n_phase)'; % time evolved during loading periods
bc.max_time_relax = max(bc.max_time_load) * ones(bc.n_phase,1); % time evolved during relaxation periods
bc.max_time = zip(bc.max_time_load',bc.max_time_relax')';
% the displacement max_displacement is applied during the time max_time
bc.max_displacement_load = 0.5*ones(bc.n_phase,1); % displacement applied during loading periods
bc.max_displacement_relax = 0.0*ones(bc.n_phase,1); % displacement applied during relaxation periods
bc.max_displacement = zip(bc.max_displacement_load',bc.max_displacement_relax')';
% time discretization
bc.n_step_load = 20*ones(bc.n_phase,1); % number of steps used during loading periods
bc.n_step_relax = 20*ones(bc.n_phase,1); % number of steps used during relaxation periods
bc.n_step = zip(bc.n_step_load',bc.n_step_relax')'; % number of load steps for each case

bc.displacement_increment = zeros(length(bc.max_displacement), length(bc.dof_displacement));
for i = 1:length(bc.max_displacement)
    if i <= length(bc.max_displacement)/2
    	bc.displacement_increment(i,:) = [zeros(size(nodes_displacement)), ones(size(nodes_displacement))]; % vertical displacement, tension
    else
        bc.displacement_increment(i,:) = [zeros(size(nodes_displacement)), -ones(size(nodes_displacement))]; % vertical displacement, compression
    end
end


