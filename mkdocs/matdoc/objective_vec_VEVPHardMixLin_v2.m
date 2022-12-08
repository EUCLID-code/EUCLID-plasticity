function cost = objective_vec_VEVPHardMixLin_v2(theta,n_Maxwell,q_case,q,lambda_r,lambda_p,p,n_eval,time,u,reaction,n_dof,n_element,n_dof_per_element,element_dof,Gauss_weights,n_Gauss_per_dim,detJ_GP,B_GP,dof_free,dof_reaction)
% objective_vec_VEVPHardMixLin_v2 calculates the objective minimization function
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
% `theta` (_double_) - material parameters
%
% `n_Maxwell` (_double_) - number of Maxwell elements
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
% `time` (_double_) - time step data
%
% `u` (_double_) - displacement data for every load step (every column
% corresponds to one load step)
%
% `reaction` (_double_) - reaction force data for every load step
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
if length(theta) == 2
    % linear elastic
    n_Maxwell_G = 0;
    n_Maxwell_K = 0;
    G0 = theta(1); K0 = theta(2);
    Gi = zeros(1,n_Maxwell_G); % Gi = []; produces a conflict in .mex file generation
    gi = zeros(1,n_Maxwell_G); % gi = []; produces a conflict in .mex file generation
    Ki = zeros(1,n_Maxwell_K); % Ki = []; produces a conflict in .mex file generation
    ki = zeros(1,n_Maxwell_K); % ki = []; produces a conflict in .mex file generation
    H_iso = 0; H_kin = 0; eta = 1e6; sigma_0 = 1e6;
elseif length(theta) == (2 + 2*sum(n_Maxwell))
    % linear viscoelastic
    n_Maxwell_G = n_Maxwell(1);
    n_Maxwell_K = n_Maxwell(2);
    counter = 1;
    G0 = theta(counter); counter = counter + 1;
    K0 = theta(counter); counter = counter + 1;
    Gi = zeros(1,n_Maxwell_G);
    gi = zeros(1,n_Maxwell_G);
    for i = 1:n_Maxwell_G
        Gi(i) = theta(counter); counter = counter + 1;
        gi(i) = 1 / theta(counter); counter = counter + 1;
    end
    Ki = zeros(1,n_Maxwell_K);
    ki = zeros(1,n_Maxwell_K);
    for i = 1:n_Maxwell_K
        Ki(i) = theta(counter); counter = counter + 1;
        ki(i) = 1 / theta(counter); counter = counter + 1;
    end
    H_iso = 0; H_kin = 0; eta = 1e6; sigma_0 = 1e6;
elseif length(theta) == (2 + 2*sum(n_Maxwell) + 4)
    % viscoelastic - viscoplastic
    % linear isotropic hardening
    % linear kinematic hardening
    % linear viscoplastic evolution
    n_Maxwell_G = n_Maxwell(1);
    n_Maxwell_K = n_Maxwell(2);
    counter = 1;
    G0 = theta(counter); counter = counter + 1;
    K0 = theta(counter); counter = counter + 1;
    Gi = zeros(1,n_Maxwell_G);
    gi = zeros(1,n_Maxwell_G);
    for i = 1:n_Maxwell_G
        Gi(i) = theta(counter); counter = counter + 1;
        gi(i) = 1 / theta(counter); counter = counter + 1;
    end
    Ki = zeros(1,n_Maxwell_K);
    ki = zeros(1,n_Maxwell_K);
    for i = 1:n_Maxwell_K
        Ki(i) = theta(counter); counter = counter + 1;
        ki(i) = 1 / theta(counter); counter = counter + 1;
    end
    H_iso = theta(counter); counter = counter + 1;
    H_kin = theta(counter); counter = counter + 1;
    eta = theta(counter); counter = counter + 1;
    sigma_0 = 1 / theta(counter);
else
    cost(:) = NaN; 
    disp('Wrong parameter input.')
    return
end
Ginf = G0 - sum(Gi);
Kinf = K0 - sum(Ki);

%% History Variables
n_component_3D = 6;
alphaV_GP_prev = zeros(n_element,n_Gauss_per_dim,n_Gauss_per_dim,n_component_3D,n_Maxwell_K);
alphaV_GP = zeros(n_element,n_Gauss_per_dim,n_Gauss_per_dim,n_component_3D,n_Maxwell_K);
epsilonVvp_GP_prev = zeros(n_element,n_Gauss_per_dim,n_Gauss_per_dim,n_component_3D);
epsilonVvp_GP = zeros(n_element,n_Gauss_per_dim,n_Gauss_per_dim,n_component_3D);
gamma_GP_prev = zeros(n_element,n_Gauss_per_dim,n_Gauss_per_dim);
gamma_GP = zeros(n_element,n_Gauss_per_dim,n_Gauss_per_dim);

%% Internal and External Forces
for idx_displacement = 1:n_eval:size(reaction,2)

    f_int_global = zeros(n_dof,1);

    if idx_displacement == 1
        time_inc = time(idx_displacement);
    else
        time_inc = time(idx_displacement) - time(idx_displacement-1);
    end

    for idx_ele = 1:n_element

        f_int_ele = zeros(n_dof_per_element,1);
        dof_ele = element_dof(idx_ele,:);
        u_ele = u(dof_ele,idx_displacement);

        for idx_Gauss_x = 1:n_Gauss_per_dim
            for idx_Gauss_y = 1:n_Gauss_per_dim

                detJ = detJ_GP(idx_ele,idx_Gauss_x,idx_Gauss_y);
                B = squeeze(B_GP(idx_ele,idx_Gauss_x,idx_Gauss_y,:,:));
                epsilonV = B*u_ele;

                % history variables at the previous time
                alphaV_prev = squeeze(alphaV_GP_prev(idx_ele,idx_Gauss_x,idx_Gauss_y,:,:));
                epsilonVvp_prev = squeeze(epsilonVvp_GP_prev(idx_ele,idx_Gauss_x,idx_Gauss_y,:));
                gamma_prev = squeeze(gamma_GP_prev(idx_ele,idx_Gauss_x,idx_Gauss_y));

                % local problem
                [sigmaV,alphaV,epsilonVvp,gamma,viscoelastic,converged_local] = ...
                    local_problem_PlaneStrain_VEVPHardMixLin_notan( ...
                    Gi, ...
                    gi, ...
                    Ki, ...
                    ki, ...
                    Ginf, ...
                    Kinf, ...
                    H_iso, ... % isotropic hardening parameter
                    H_kin, ... % kinematic hardening parameter
                    eta, ...
                    sigma_0, ...
                    time_inc, ...
                    epsilonV, ... % strain increment
                    alphaV_prev, ... % viscoelastic strain @ n
                    epsilonVvp_prev, ... % viscoplastic strain @ n
                    gamma_prev ... % plastic multiplier @ n
                    );
                if ~converged_local
                    cost(:) = NaN; %disp('cost = NaN')
                    return
                end

                % history variables at the current time
                alphaV_GP(idx_ele,idx_Gauss_x,idx_Gauss_y,:,:) = alphaV;
                epsilonVvp_GP(idx_ele,idx_Gauss_x,idx_Gauss_y,:) = epsilonVvp;
                gamma_GP(idx_ele,idx_Gauss_x,idx_Gauss_y) = gamma;

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

    alphaV_GP_prev = alphaV_GP;
    epsilonVvp_GP_prev = epsilonVvp_GP;
    gamma_GP_prev = gamma_GP;

end

%% Regularization
cost(end) = sqrt(lambda_p*sum(abs(theta(3:end)).^p));

end










