function cost = objective_vec_HardMixVAF(parameters,q_case,q,lambda_r,lambda_p,p,n_eval,u,reaction,n_NR_local,tol_NR_local,CPlaneStress,SPlaneStress,n_dof,n_element,n_dof_per_element,element_dof,Gauss_weights,n_Gauss_per_dim,detJ_GP,B_GP,dof_free,dof_reaction)
% objective_vec_HardMixVAF calculates the objective minimization function
%
% ## Comments 
% 
% This function defines the objective function that depends on
% the given data and on theta. Note that the function
% does not return a scalar but a vector. The Matlab build-in optimizer
% (e.g., lsqnonlin) minimizes the sum of squared entries of this vector.
%
% This function may be converted into a __Matlab executable file__ (.mex
% file) with the [__Matlab Coder__](https://mathworks.com/help/coder/) to
% accelerate the objective function evaluation.
%
% ## Input Arguments
% 
% `paramters` (_double_) - material parameters
% 
% `q_case` (_double_) - see `opt_input` documentation
% 
% `q` (_double_) - see `opt_input` documentation
% 
% `lambda_r` (_double_) - see `opt_input` documentation
% 
% `lambda_p` (_double_) - see `opt_input` documentation
% 
% `p` (_double_) - see `opt_input` documentation
% 
% `n_eval` (_double_) - see `opt_input` documentation
% 
% `u` (_double_) - displacement data for every load step (every column
% corresponds to one load step)
% 
% `reaction` (_double_) - reaction force data for every load step
% 
% `n_NR_local` (_double_) - maximum number of Newton-Raphson iterations
% used for solving the local problem (return mapping algorithm)
% 
% `tol_NR_local` (_double_) - stopping tolerance for the Newton-Raphson
% iteration used for solving the local problem (return mapping algorithm)
% 
% `CPlaneStress` (_double_) - known elastic properties (stiffness matrix)
% 
% `SPlaneStress` (_double_) - known elastic properties (compliance matrix)
% 
% `n_dof` (_double_) - number of degrees of freedom in the mesh
% 
% `n_element` (_double_) - number of finite elements in the mesh
% 
% `n_dof_per_element` (_double_) - number of degrees of freedom per finite
% element
% 
% `element_dof` (_double_) - degrees of freedom at each finite element
% element (each row provides all degrees of freedom at the corresponding
% finite element)
% 
% `Gauss_weights` (_double_) - Gaussian quadrature weights
% 
% `n_Gauss_per_dim` (_double_) - number of Gauss points per dimension
% 
% `detJ_GP` (_double_) - determinant of the Jacobian at each Gauss point
% needed for mapping integrals between the reference and physical finite
% element
% 
% `B_GP` (_double_) - help matrix at each Gauss point needed for
% calculating shape function derivatives
% 
% `dof_free` (_double_) - free degrees of freedom
% 
% `dof_reaction` (_double_) - degrees of freedom corresponding to reaction
% forces
% 
% ## Output Arguments
% 
% `cost` (_double_) - vector, whose sum of squared entries must be
% minimized (note that the regularization term is included as the last
% entry of this vector)

%% Allocation
n_eval_step = length(1:n_eval:size(reaction,2));
n_residual = length(dof_free) + size(dof_reaction,1);
cost = zeros(n_eval_step*n_residual+1,1);

%% Material Parameters
if length(parameters) == 1
    theta = parameters;
    H_isotropic = zeros(1,3);
    H_kinematic = zeros(1,2);
elseif length(parameters) == 3
    theta = parameters(1);
    H_isotropic = [parameters(2), 0, 0];
    H_kinematic = [parameters(3), 0];
elseif length(parameters) > 5
    theta = parameters(1:end-5);
    H_isotropic = parameters(end-4:end-2);
    H_kinematic = parameters(end-1:end);
else
    cost(:) = NaN; 
    disp('Wrong parameter input.') %disp('cost = NaN')
    return
end


%% Physical Requirements
if theta(1) <= sum(abs(theta(2:end)))
    cost(:) = NaN;
%     disp('Fourier coefficients are not admissible.') %disp('cost = NaN')
    return
end

%% Hystory Variables
epsilonVp_GP_prev = zeros(n_element,n_Gauss_per_dim,n_Gauss_per_dim,4);
epsilonVp_GP = zeros(n_element,n_Gauss_per_dim,n_Gauss_per_dim,4);
gamma_GP_prev = zeros(n_element,n_Gauss_per_dim,n_Gauss_per_dim);
gamma_GP = zeros(n_element,n_Gauss_per_dim,n_Gauss_per_dim);
sigmaV_back_GP_prev = zeros(n_element,n_Gauss_per_dim,n_Gauss_per_dim,3);
sigmaV_back_GP = zeros(n_element,n_Gauss_per_dim,n_Gauss_per_dim,3);

%% Internal and External Forces
for idx_displacement = 1:n_eval:size(reaction,2)
        
    f_int_global = zeros(n_dof,1);
    
    for idx_ele = 1:n_element
        
        f_int_ele = zeros(n_dof_per_element,1);
        dof_ele = element_dof(idx_ele,:);
        u_ele = u(dof_ele,idx_displacement);
        
        for idx_Gauss_x = 1:n_Gauss_per_dim
            for idx_Gauss_y = 1:n_Gauss_per_dim
                
                detJ = detJ_GP(idx_ele,idx_Gauss_x,idx_Gauss_y);
                B = squeeze(B_GP(idx_ele,idx_Gauss_x,idx_Gauss_y,:,:));
                epsilonV = B*u_ele;
                epsilonVp_prev = squeeze(epsilonVp_GP_prev(idx_ele,idx_Gauss_x,idx_Gauss_y,:));
                gamma_prev = gamma_GP_prev(idx_ele,idx_Gauss_x,idx_Gauss_y);
                sigmaV_back_prev = squeeze(sigmaV_back_GP_prev(idx_ele,idx_Gauss_x,idx_Gauss_y,:));
                
                [sigmaV,epsilonVp,gamma,sigmaV_back,converged_local] = local_problem_PlaneStressHardMixVAF_notan( ...
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
                if ~converged_local
                    cost(:) = NaN; %disp('cost = NaN')
                    return
                end
                epsilonVp_GP(idx_ele,idx_Gauss_x,idx_Gauss_y,:) = epsilonVp;
                gamma_GP(idx_ele,idx_Gauss_x,idx_Gauss_y) = gamma;
                sigmaV_back_GP(idx_ele,idx_Gauss_x,idx_Gauss_y,:) = sigmaV_back;
                
                f_int_ele = f_int_ele ...
                    + B' * sigmaV * Gauss_weights(idx_Gauss_x) * Gauss_weights(idx_Gauss_y) * detJ;
                
            end
        end
        
        % assembly
        f_int_global(dof_ele) = f_int_global(dof_ele) + f_int_ele;
        
    end
    
    residual_free = f_int_global(dof_free);
    residual_fix = sum(f_int_global(dof_reaction),2) - reaction(:,idx_displacement);
        
    cost_idx_displacment = [residual_free; sqrt(lambda_r)*residual_fix];
    switch q_case
        case 0
            if q == 2
                
            else
                disp('Not yet implemented.')
            end
        case 1
            disp('Not yet implemented.')
        case 2
            disp('Not yet implemented.')
    end
    
    cost(  ((idx_displacement-1)*n_residual+1)  :  (idx_displacement*n_residual)  ) = cost_idx_displacment;
    
    epsilonVp_GP_prev = epsilonVp_GP;
    gamma_GP_prev = gamma_GP;
    sigmaV_back_GP_prev = sigmaV_back_GP;
    
end

%% Regularization
if length(theta) > 1
    cost(end) = sqrt(lambda_p*sum(abs(theta(2:end)).^p));
else
    cost(end) = 0;
end

end










