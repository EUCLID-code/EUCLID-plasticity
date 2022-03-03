function FEA_plate_elliptic_holes()
% FEA_plate_elliptic_holes executes finite element analysis
%
% ## Comments
% 
% Finite element data is saved as a structural array. The material, mesh
% and boundary conditions have to be specified by the user. Set `testcase = true` and `material = vonMises` for a computationally inexpensive example.
%
% Assumptions:
% - two-dimensional
% - plane stress
% - no body force
% - only homogeneous Neumann boundary conditions
% - displacement control
% 
% ## Input Arguments
% 
% _none_
% 
% ## Output Arguments
% 
% _none_
% 

addpath(genpath(pwd));
testcase = true;

%% Algorithm
% Newton-Raphson
algorithm.n_NR_global = 100; % maximum number of global Newton-Raphson iterations
algorithm.tol_NR_global = 1e-9; % convergence tolerance global Newton-Raphson algorithm
algorithm.n_NR_local = 100; % maximum number of local Newton-Raphson iterations
algorithm.tol_NR_local = 1e-10; % convergence tolerance local Newton-Raphson algorithm
% save stress states at which plastic deformation occurs, i.e., at which
% the yield stress is reached
algorithm.n_save_YieldSurface = 3600;

%% Material
material = vonMises;
disp(['Material: ' material.name])

%% Mesh
% (i) if the mesh has not been used before:
% mesh.name = 'plate_elliptic_holes';
% L = 100;
% load('node_plate_elliptic_holes_x1.mat'); node = L*node; mesh.node = node;
% load('element_plate_elliptic_holes_x1.mat'); mesh.element = element;
% mesh.bot_edge = find(abs(node(:,2)) < 1e-8)';
% mesh.top_edge = find(abs(node(:,2) - L) < 1e-8)';
% mesh = mesh_info(mesh); save('mesh/mesh_plate_elliptic_holes_x1.mat','mesh')
% (ii) else:
load('mesh_plate_elliptic_holes_x1.mat','mesh'); disp(['Mesh: ' mesh.name])

%% Boundary Conditions
if testcase
    bc = bc_plate_elliptic_holes_UT_testcase(mesh);
else
    bc = bc_plate_elliptic_holes_UTC5(mesh);
%     bc = bc_plate_elliptic_holes_shear(mesh);
end
disp(['Boundary condition: ' bc.name])

%% Global Problem (inherits Local Problem)
results = global_problem_call(algorithm,material,mesh,bc);

%% Save Data
save_path = [datestr(now,'yyyymmddTHHMMSS') '_DATA_' material.name '_' mesh.name];
save_path = ['FEM_data\' save_path];
data.algorithm = algorithm; data.material = material; data.mesh = mesh; data.bc = bc; data.results = results;
save(save_path,'data')

end
