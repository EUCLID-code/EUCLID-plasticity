function [time,u,reaction,alphaV_GP,epsilonVvp_GP,gamma_GP] = global_problem_PlaneStrain_VEVPHardMixLin( ...
    ... % algorithm
    n_NR_global, ...
    tol_NR_global, ...
    ... % material
    material, ...
    ... % mesh
    n_dof, ...
    n_element, ...
    n_dof_per_element, ...
    element_dof, ...
    Gauss_weights, ...
    n_Gauss_per_dim, ...
    detJ_GP, ...
    B_GP, ...
    ... % bc
    dof_fix, ...
    dof_displacement, ...
    dof_reaction, ...
    max_time, ...
    max_displacement, ...
    n_step, ...
    displacement_increment ...
    )

% ==================================================
% Solving the global problem.
% ==================================================

%% Initial Condition
u_global = zeros(n_dof,1);
total_n_load = sum(n_step);
displacement = zeros(1,size(displacement_increment,2));

%% Load
n_reaction = size(dof_reaction,1);
n_load_case = length(max_displacement);

%% History Variables
n_component_3D = 6;
if material.n_Maxwell_G ~= material.n_Maxwell_K
    error('Not yet implemented.');
end
alphaV_GP_prev = zeros(n_element,n_Gauss_per_dim,n_Gauss_per_dim,n_component_3D,material.n_Maxwell_K);
alphaV_GP = zeros(n_element,n_Gauss_per_dim,n_Gauss_per_dim,n_component_3D,material.n_Maxwell_K);
epsilonVvp_GP_prev = zeros(n_element,n_Gauss_per_dim,n_Gauss_per_dim,n_component_3D);
epsilonVvp_GP = zeros(n_element,n_Gauss_per_dim,n_Gauss_per_dim,n_component_3D);
gamma_GP_prev = zeros(n_element,n_Gauss_per_dim,n_Gauss_per_dim);
gamma_GP = zeros(n_element,n_Gauss_per_dim,n_Gauss_per_dim);

%% Allocation
counter = 1;
% f_int_global = zeros(n_dof,1);
% K_global = zeros(n_dof);
residual_norm_0 = 1;
time = zeros(1,total_n_load);
u = zeros(n_dof,total_n_load);
reaction = zeros(n_reaction,total_n_load);

%% Solve
for idx_displacement_case = 1:n_load_case
    
    fprintf('\n')
    fprintf('Displacement Case: %i\n',idx_displacement_case)
    
    time_inc = max_time(idx_displacement_case)/n_step(idx_displacement_case);

    for idx_displacement = 1:n_step(idx_displacement_case)
        
        fprintf('Displacement Step: %i\n',idx_displacement)
        
        displacement = displacement + displacement_increment(idx_displacement_case,:)*max_displacement(idx_displacement_case)/n_step(idx_displacement_case);
        
        dof_Dirichlet = [dof_fix, dof_displacement];
        value_Dirichlet = [zeros(size(dof_fix)), displacement];
        
        converged_global = false;
        
        for idx_NR = 1:n_NR_global
            
            f_int_global = zeros(n_dof,1);
            K_global = zeros(n_dof);
            
            for idx_ele = 1:n_element
                
                f_int_ele = zeros(n_dof_per_element,1);
                K_ele = zeros(n_dof_per_element);
                dof_ele = element_dof(idx_ele,:);
                u_ele = u_global(dof_ele);
                
                for idx_Gauss_x = 1:n_Gauss_per_dim
                    for idx_Gauss_y = 1:n_Gauss_per_dim
                        
                        % strain increment
                        detJ = detJ_GP(idx_ele,idx_Gauss_x,idx_Gauss_y);
                        B = squeeze(B_GP(idx_ele,idx_Gauss_x,idx_Gauss_y,:,:));
                        epsilonV = B*u_ele;

                        % history variables at the previous time
                        alphaV_prev = squeeze(alphaV_GP_prev(idx_ele,idx_Gauss_x,idx_Gauss_y,:,:));
                        epsilonVvp_prev = squeeze(epsilonVvp_GP_prev(idx_ele,idx_Gauss_x,idx_Gauss_y,:));
                        gamma_prev = squeeze(gamma_GP_prev(idx_ele,idx_Gauss_x,idx_Gauss_y));

                        % local problem
                        [sigmaV,CC,alphaV,epsilonVvp,gamma,viscoelastic,converged_local] = local_problem_PlaneStrain_VEVPHardMixLin( ...
                            material.Gi, ...
                            material.gi, ...
                            material.Ki, ...
                            material.ki, ...
                            material.Ginf, ...
                            material.Kinf, ...
                            material.H_iso, ... % isotropic hardening parameter
                            material.H_kin, ... % kinematic hardening parameter
                            material.eta, ...
                            material.sigma_0, ...
                            time_inc, ...
                            epsilonV, ... % strain increment
                            alphaV_prev, ... % viscoelastic strain @ n
                            epsilonVvp_prev, ... % viscoplastic strain @ n
                            gamma_prev ... % plastic multiplier @ n
                            );
                        if ~converged_local
                            return
                        end

                        % history variables at the current time
                        alphaV_GP(idx_ele,idx_Gauss_x,idx_Gauss_y,:,:) = alphaV;
                        epsilonVvp_GP(idx_ele,idx_Gauss_x,idx_Gauss_y,:) = epsilonVvp;
                        gamma_GP(idx_ele,idx_Gauss_x,idx_Gauss_y) = gamma;
                        
                        % internal forces and stiffness matrices at the
                        % element
                        f_int_ele = f_int_ele ...
                            + B' * sigmaV * Gauss_weights(idx_Gauss_x) * Gauss_weights(idx_Gauss_y) * detJ;
                        K_ele = K_ele + ...
                            B' * CC * B * Gauss_weights(idx_Gauss_x) * Gauss_weights(idx_Gauss_y) * detJ;
                        
                    end
                end
                
                % assembly
                f_int_global(dof_ele) = f_int_global(dof_ele) + f_int_ele;
                K_global(dof_ele,dof_ele) = K_global(dof_ele,dof_ele) + K_ele;
                
            end
            
            % reaction forces
            reaction(:,counter) = sum(f_int_global(dof_reaction),2);
            
            % Dirichlet boundary conditions
            for d = 1:length(dof_Dirichlet)
                K_global(dof_Dirichlet(d),:) = 0;
                K_global(dof_Dirichlet(d),dof_Dirichlet(d)) = 1;
                % no external force (f_ext_global)
                f_int_global(dof_Dirichlet(d)) = ...
                    u_global(dof_Dirichlet(d)) - value_Dirichlet(d);
            end
            
            % convergence
            residual_norm = norm(f_int_global);
            if idx_NR == 1
                residual_norm_0 = residual_norm;
                fprintf('Order of residual: ')
%                 fprintf('Order of relative residual: ')
            end
            fprintf('%d, ',ceil(log10(residual_norm)))
%             fprintf('%d ',ceil(log10(residual_norm/residual_norm_0)))
            if (residual_norm / residual_norm_0) < tol_NR_global || residual_norm < 1e-2*tol_NR_global
                fprintf('convergence after %d iterations.\n',idx_NR-1)
                converged_global = true;
                break
            end
            
            % update displacements
            u_global = u_global - K_global\f_int_global;
            
        end
        
        if converged_global
            alphaV_GP_prev = alphaV_GP;
            epsilonVvp_GP_prev = epsilonVvp_GP;
            gamma_GP_prev = gamma_GP;
            if counter == 1
                time(1,counter) = time_inc;
            else
                time(1,counter) = time(1,counter-1) + time_inc;
            end
            u(:,counter) = u_global;
            counter = counter + 1;
        else
            fprintf('Newton-Raphson of the global problem did not converge.\n')
            return
        end
        
    end
    
end

end

