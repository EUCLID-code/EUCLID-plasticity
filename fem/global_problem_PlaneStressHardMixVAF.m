function [u,reaction,epsilonVp_GP,gamma_GP,sigmaV_back_GP,YieldSurface] = global_problem_PlaneStressHardMixVAF(n_save_YieldSurface,n_NR_global,tol_NR_global,n_NR_local,tol_NR_local,CPlaneStress,SPlaneStress,theta,H_isotropic,H_kinematic,n_dof,n_element,n_dof_per_element,element_dof,Gauss_weights,n_Gauss_per_dim,detJ_GP,B_GP,dof_fix,dof_displacement,dof_reaction,max_displacement,n_step,displacement_increment)
% global_problem_PlaneStressHardMixVAF solves the global finite element problem
%
% ## Comments
% 
% This function solves the global finite element problem, i.e., the
% residuals of the discretized weak formulation are minimized. Assumptions:
% plane stress, isotropic hardening (Voce), kinematic hardening
% (Armstrong-Frederick).
% 
% ## Input Arguments
% 
% `n_save_YieldSurface` (_double_) - 
% 
% `n_NR_global` (_double_) - maximum number of Newton-Raphson iterations
% used for solving the global problem
% 
% `tol_NR_global` (_double_) - stopping tolerance for the Newton-Raphson
% iteration used for solving the global problem
% 
% `n_NR_local` (_double_) - maximum number of Newton-Raphson iterations
% used for solving the local problem (return mapping algorithm)
% 
% `tol_NR_local` (_double_) - stopping tolerance for the Newton-Raphson
% iteration used for solving the local problem (return mapping algorithm)
% 
% `CPlaneStress` (_double_) - elastic properties (stiffness matrix)
% 
% `SPlaneStress` (_double_) - elastic properties (compliance matrix)
% 
% `theta` (_double_) - material parameters
%
% `H_isotropic` (_double_) - isotropic hardening parameters
%
% `H_kinematic` (_double_) - kinematic hardening parameters
% 
% `n_dof` (_double_) - see `mesh_info`
% 
% `n_element` (_double_) - see `mesh_info`
% 
% `n_dof_per_element` (_double_) - see `mesh_info`
% 
% `element_dof` (_double_) - see `mesh_info`
% 
% `Gauss_weights` (_double_) - see `mesh_info`
% 
% `n_Gauss_per_dim` (_double_) - see `mesh_info`
% 
% `detJ_GP` (_double_) - see `mesh_info`
% 
% `B_GP` (_double_) - see `mesh_info`
% 
% `dof_fix` (_double_) - fixed degrees of freedom
% 
% `dof_displacement` (_double_) - degrees of freedom under displacement
% control
% 
% `dof_reaction` (_double_) - degrees of freedom associated with the
% reaction forces
% 
% `max_displacement` (_double_) - maximum applied displacement
% 
% `n_step` (_double_) - number of load steps
% 
% `displacement_increment` (_double_) - applied displacement increment
% 
% ## Output Arguments
% 
% `u` (_double_) - displacement at each node for each load step
%
% `reaction` (_double_) - reaction forces for each load step
%
% `epsilonVp_GP` (_double_) - plastic component of the infinitesimal
% strain at each Gauss point for the last load step
%
% `gamma_GP` (_double_) - plastic multiplier at each Gauss point for the last load step
%
% `sigmaV_back_GP` (_double_) - back stress at each Gauss point for the last load step
%
% `YieldSurface` (_double_) - stresses at which yielding has occured 
%

u_global = zeros(n_dof,1); % initial condition
total_n_load = sum(n_step);
displacement = zeros(1,size(displacement_increment,2));

%% Load
n_reaction = size(dof_reaction,1);
n_load_case = length(max_displacement);

%% History Variables
epsilonVp_GP_prev = zeros(n_element,n_Gauss_per_dim,n_Gauss_per_dim,4);
epsilonVp_GP = zeros(n_element,n_Gauss_per_dim,n_Gauss_per_dim,4);
gamma_GP_prev = zeros(n_element,n_Gauss_per_dim,n_Gauss_per_dim);
gamma_GP = zeros(n_element,n_Gauss_per_dim,n_Gauss_per_dim);
sigmaV_back_GP_prev = zeros(n_element,n_Gauss_per_dim,n_Gauss_per_dim,3);
sigmaV_back_GP = zeros(n_element,n_Gauss_per_dim,n_Gauss_per_dim,3);

%% Allocation
counter = 1;
% f_int_global = zeros(n_dof,1);
% K_global = zeros(n_dof);
residual_norm_0 = 1;
u = zeros(n_dof,total_n_load);
reaction = zeros(n_reaction,total_n_load);
YieldSurface = zeros(n_save_YieldSurface,2);

%% Solve
for idx_displacement_case = 1:n_load_case
    
    %     fprintf('Displacement Case: %i\n',idx_displacement_case)
    
    for idx_displacement = 1:n_step(idx_displacement_case)
        
        %         fprintf('Displacement Step: %i\n',idx_displacement)
        
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
                        
                        detJ = detJ_GP(idx_ele,idx_Gauss_x,idx_Gauss_y);
                        B = squeeze(B_GP(idx_ele,idx_Gauss_x,idx_Gauss_y,:,:));
                        epsilonV = B*u_ele;
                        epsilonVp_prev = squeeze(epsilonVp_GP_prev(idx_ele,idx_Gauss_x,idx_Gauss_y,:));
                        gamma_prev = gamma_GP_prev(idx_ele,idx_Gauss_x,idx_Gauss_y);
                        sigmaV_back_prev = squeeze(sigmaV_back_GP_prev(idx_ele,idx_Gauss_x,idx_Gauss_y,:));
                        
                        [sigmaV,CC,epsilonVp,elastic,gamma,sigmaV_back,converged_local] = local_problem_PlaneStressHardMixVAF( ...
                            n_NR_local, ...
                            tol_NR_local, ...
                            CPlaneStress, ...
                            SPlaneStress, ...
                            theta, ...
                            H_isotropic, ...
                            H_kinematic, ...
                            epsilonV, ...
                            epsilonVp_prev, ...
                            gamma_prev, ...
                            sigmaV_back_prev ...
                            );
                        
                        if ~elastic
                            [r, alpha] = get_r_alpha_PlaneStress(sigmaV - sigmaV_back);
                            position = floor(alpha*n_save_YieldSurface/(2*pi));
                            while position < 0
                                position = position + n_save_YieldSurface;
                            end
                            while position > n_save_YieldSurface
                                position = position - n_save_YieldSurface;
                            end
                            YieldSurface(position,:) = [alpha, r];
                        end
                        if ~converged_local
                            return
                        end
                        epsilonVp_GP(idx_ele,idx_Gauss_x,idx_Gauss_y,:) = epsilonVp;
                        gamma_GP(idx_ele,idx_Gauss_x,idx_Gauss_y) = gamma;
                        sigmaV_back_GP(idx_ele,idx_Gauss_x,idx_Gauss_y,:) = sigmaV_back;
                        
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
%                 fprintf('Order of Residual: ')
            end
%             fprintf('%d ',ceil(log10(residual_norm/residual_norm_0)))
            if (residual_norm / residual_norm_0) < tol_NR_global
%                 fprintf('\nConvergence after %d iterations.\n',idx_NR-1)
                converged_global = true;
                break
            end
            
            % update displacements
            u_global = u_global - K_global\f_int_global;
            
        end
        
        if converged_global
            epsilonVp_GP_prev = epsilonVp_GP;
            gamma_GP_prev = gamma_GP;
            sigmaV_back_GP_prev = sigmaV_back_GP;
            u(:,counter) = u_global;
            counter = counter + 1;
        else
            disp('Newton-Raphson of the global problem did not converge.')
            break
        end
        
    end
    
end

end

