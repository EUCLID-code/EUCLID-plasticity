function opt_HardMix()
% opt_HardMix runs the optimization problem
%
% ## Comments
%
% Please make sure that the correct input is defined in
% `opt_input_HardMix`. The results are saved in the directory specified in
% `opt_input_HardMix`.
%
% ## Input Arguments
%
% _none_
%
% ## Output Arguments
%
% _none_

in = opt_input_HardMix(); % load input paramaters / hyperparamters / settings etc.

%% Load Data
load(in.datafile,'data')
disp(['Material: ' data.material.name])
disp(['Mesh: ' data.mesh.name])
disp(['Boundary Condition: ' data.bc.name])
disp(['Data: ' in.datafile])

%% Optimization Method
options = optimset( ...
    'MaxFunEvals',in.MaxFunEvals, ...
    'TolFun',in.TolFun, ...
    'TolX',in.TolX, ...
    'Display','off' ... % 'off' or 'iter'
    );

tic
disp(['Method: ' in.method])
% to solve the optimization problem, different optimization algorithms may
% be used
% here, only one algorithm 'lsqnonlin_lb_parsplit_HardMixVAF' is
% implemented
if strcmp(in.method,'lsqnonlin_lb_parsplit_HardMixVAF')
    %% Objective Function without Regularization
    in_without_reg = in;
    in_without_reg.lambda_p = 0;
    % define objective function without regularization term
    fun_vec = def_objective_vec_HardMixVAF(in_without_reg,data);

    %% Test True Material Parameters
    cost_true = Inf;
    cost_true_crop = Inf;
    if contains(data.material.name,'FYS')
        % it is interesting to compute the
        % objective function for the true parameter vector
        theta_true = data.material.theta;
        H_isotropic_true = data.material.H_isotropic;
        H_kinematic_true = data.material.H_kinematic;
        parameters_true = [theta_true, H_isotropic_true, H_kinematic_true];
        disp('Testing true material parameters:')
        disp(['Parameters (theta): ' num2str(theta_true)])
        disp(['Parameters (H_isotropic): ' num2str(H_isotropic_true)])
        disp(['Parameters (H_kinematic): ' num2str(H_kinematic_true)])
        tic
        % compute the objective function for the true parameters
        cost_true = fun_vec(parameters_true);
        cost_true = sum(cost_true.^2);
        disp(['Cost: ' num2str(cost_true)])
        toc % measure time to compute the objective function

        if length(data.material.theta) > in.n_feature
            % if the true material parameter vector theta has more elements
            % than the number of features chosen, the true solution cannot
            % be element of the search space, i.e., it is not possible to
            % calculate the exact material parameters (even for the case
            % without noise)
            % in this case, it is interesting to compute the
            % objective function for the true parameter vector, which is
            % cropped such that the number of parameters equals the number
            % of features chosen
            theta_true_crop = data.material.theta(1:in.n_feature);
            parameters_true_crop = [theta_true_crop, H_isotropic_true, H_kinematic_true];
            disp('Testing true material parameters:')
            disp(['Parameters (theta): ' num2str(theta_true_crop)])
            disp(['Parameters (H_isotropic): ' num2str(H_isotropic_true)])
            disp(['Parameters (H_kinematic): ' num2str(H_kinematic_true)])
            tic
            % compute the objective function for the cropped true parameter
            % vector
            cost_true_crop = fun_vec(parameters_true_crop);
            cost_true_crop = sum(cost_true_crop.^2);
            disp(['Cost: ' num2str(cost_true_crop)])
            toc % measure time to compute the objective function
        end
    end

    %% Solution Array
    solution.in = in;
    if contains(data.material.name,'FYS')
        solution.theta_true = theta_true;
        solution.H_isotropic_true = H_isotropic_true;
        solution.H_kinematic_true = H_kinematic_true;
        solution.cost_true = cost_true;
        if length(data.material.theta) > in.n_feature
            solution.cost_true_crop = cost_true_crop;
        end
    end
    if ~contains(in.datafile,'noise')
        solution.in.standard_deviation = 0;
    else
        solution.in.standard_deviation = data.results.standard_deviation;
    end

    %% Initial Guess
    disp('Initial Guess:')
    theta_init = in.theta_init;
    theta_save = theta_init;
    H_isotropic_init = in.H_isotropic_init;
    H_isotropic_save = H_isotropic_init;
    H_kinematic_init = in.H_kinematic_init;
    H_kinematic_save = H_kinematic_init;
    parameters_init = [theta_init, H_isotropic_init, H_kinematic_init];
    disp(['Solution (theta): ' num2str(theta_init)])
    disp(['Solution (H_isotropic): ' num2str(H_isotropic_init)])
    disp(['Solution (H_kinematic): ' num2str(H_kinematic_init)])

    tic
    % compute the objective function for the initial parameter guess
    cost_init = fun_vec(parameters_init);
    cost_init = sum(cost_init.^2);
    cost_save = cost_init;
    disp(['Cost: ' num2str(cost_init)])
    toc % measure time to compute the objective function

    %% Preconditioner
    disp(['Preconditioner: ' in.preconditioner])
    % to accelerate the optimization process, different preconditioning
    % steps may be used
    % the idea is to solve the problem with a reduced set of parameters
    % first and then use the result as the initial condition for the higher
    % dimensional optimization
    if strcmp(in.preconditioner,'vonMises')
        [theta_init(1), min_cost, ~, ~, ~] = lsqnonlin(fun_vec,theta_init(1),[],[],options);
        theta_save = [theta_save; theta_init];
        H_isotropic_save = [H_isotropic_save; [0, 0, 0]];
        H_kinematic_save = [H_kinematic_save; [0, 0]];
    elseif strcmp(in.preconditioner,'vonMises_HardMixLin') || strcmp(in.preconditioner,'vonMises_HardMixLin_theta_0')
        if contains(in.method,'_lb_')
            [parameters_init, min_cost, ~, ~, ~] = lsqnonlin(fun_vec,[theta_init(1), H_isotropic_init(1), H_kinematic_init(1)],[-Inf, 0, 0],[],options);
        else
            [parameters_init, min_cost, ~, ~, ~] = lsqnonlin(fun_vec,[theta_init(1), H_isotropic_init(1), H_kinematic_init(1)],[],[],options);
        end
        theta_init(1) = parameters_init(1);
        H_isotropic_init(1) = parameters_init(2);
        H_kinematic_init(1) = parameters_init(3);
        theta_save = [theta_save; theta_init];
        H_isotropic_save = [H_isotropic_save; [H_isotropic_init(1), 0, 0]];
        H_kinematic_save = [H_kinematic_save; [H_kinematic_init(1), 0]];
        if strcmp(in.preconditioner,'vonMises_HardMixLin_theta_0')
            % use only theta_0 as preconditioner
            H_isotropic_init(1) = 0;
            H_kinematic_init(1) = 0;
        end
    end
    disp(['Solution (theta): ' num2str(theta_save(end,:))])
    disp(['Solution (H_isotropic): ' num2str(H_isotropic_save(end,:))])
    disp(['Solution (H_kinematic): ' num2str(H_kinematic_save(end,:))])
    cost_save = [cost_save; min_cost];
    disp(['Cost: ' num2str(min_cost)])

    %% Random Initial Guesses
    % generate random initial guesses by perturbing the solution obtained
    % from the preconditioning step
    theta_random = repmat(theta_init,in.n_guess,1);
    H_isotropic_random = repmat(H_isotropic_init,in.n_guess,1);
    H_kinematic_random = repmat(H_kinematic_init,in.n_guess,1);
    for idx = 2:in.n_guess
        admissible = false;
        while ~admissible
            theta_random(idx,:) = theta_init + in.perturbation_guess_theta.*randn(1,in.n_feature);
            % check if random initial guess is admissible
            if theta_random(idx,1) > sum(abs(theta_random(idx,2:end)))
                admissible = true;
            end
        end
        % random perturbations of the hardening parameters are assumed to be positive
        H_isotropic_random(idx,:) = H_isotropic_init + in.perturbation_guess_H_isotropic.*abs(randn(1,length(H_isotropic_init)));
        H_kinematic_random(idx,:) = H_kinematic_init + in.perturbation_guess_H_kinematic.*abs(randn(1,length(H_kinematic_init)));
    end
    solution.theta_random = theta_random;
    solution.H_isotropic_random = H_isotropic_random;
    solution.H_kinematic_random = H_kinematic_random;
    parameters_random = [theta_random, H_isotropic_random, H_kinematic_random];

    %% Parallel Optimization without Regularization
    % minimize the objective function without regularization term for the
    % randomly generated initial guesses
    % optimizations for different initial guesses are executed in parallel
    parameters_random_sol = zeros(size(parameters_random));
    cost_random_sol = zeros(in.n_guess,1);
    cluster = parcluster('local');
    pool = parpool(cluster,in.n_workers);

    disp(['Total number of parallel optimizations with different initial guesses:' num2str(in.n_guess)])
    n_parameters = size(parameters_random_sol,2);
    lb = -Inf(1,n_parameters); % lower bound
    lb(end-4:end) = 0; % hardening parameters are assumed to be positive
    parfor idx = 1:in.n_guess
        try
            [parameters, min_cost, ~, ~, ~] = lsqnonlin(fun_vec,parameters_random(idx,:),lb,[],options);
            parameters_random_sol(idx,:) = parameters;
            cost_random_sol(idx) = min_cost;
            disp(['Parallel optimization number ' num2str(idx) ' converged.'])
        catch
            % if lsqnonlin fails to solve the problem, the parameters are
            % set to NaN
            parameters_random_sol(idx,:) = NaN;
            cost_random_sol(idx) = Inf;
            disp(['Parallel optimization number ' num2str(idx) ' aborted.'])
        end
    end

    %% Solution
    theta_random_sol = parameters_random_sol(:,1:end-5);
    H_isotropic_random_sol = parameters_random_sol(:,end-4:end-2);
    H_kinematic_random_sol = parameters_random_sol(:,end-1:end);
    theta_save = [theta_save; theta_random_sol];
    H_isotropic_save = [H_isotropic_save; H_isotropic_random_sol];
    H_kinematic_save = [H_kinematic_save; H_kinematic_random_sol];
    cost_save = [cost_save; cost_random_sol];

    %% Solution Selection
    % select the solution with the lowest objective function value
    % the selected solution will later be used as the initial guess for the
    % regularized optimization
    [~,pos] = min(cost_random_sol);
    theta_init = theta_random_sol(pos,:);
    H_isotropic_init = H_isotropic_random_sol(pos,:);
    H_kinematic_init = H_kinematic_random_sol(pos,:);
    parameters_init = [theta_init, H_isotropic_init, H_kinematic_init];

    %% Regularization
    lambda_p_set = zeros(in.n_lambda_p,1);
    for idx = 1:in.n_lambda_p
        lambda_p_set(idx) = in.lambda_p*in.lambda_p_factor^(idx-1);
    end
    solution.lambda_p_set = [zeros(size(cost_save)); lambda_p_set];

    %% Parallel Optimization with Regularization
    % minimize the objective function with regularization term for
    % different choices of the regularization factor lambda_p
    % optimizations for different choices of lambda_p are executed in
    % parallel
    parameters_reg_sol = zeros(in.n_lambda_p,length(parameters_init));
    cost_reg_sol = zeros(in.n_lambda_p,1);

    disp(['Total number of parallel optimizations with varying regularization:' num2str(in.n_lambda_p)])
    n_parameters = size(parameters_random_sol,2);
    lb = -Inf(1,n_parameters); % lower bound
    lb(end-4:end) = 0; % hardening parameters are assumed to be positive
    parfor idx = 1:in.n_lambda_p
        in_mod = in;
        in_mod.lambda_p = lambda_p_set(idx);
        fun_vec = def_objective_vec_HardMixVAF(in_mod,data);
        try
            [parameters, min_cost, ~, ~, ~] = lsqnonlin(fun_vec,parameters_init,lb,[],options);
            parameters_reg_sol(idx,:) = parameters;
            cost_reg_sol(idx) = min_cost;
            disp(['Parallel optimization number ' num2str(idx) ' converged.'])
        catch
            % if lsqnonlin fails to solve the problem, the parameters are
            % set to NaN
            parameters_reg_sol(idx,:) = NaN;
            cost_reg_sol(idx) = Inf;
            disp(['Parallel optimization number ' num2str(idx) ' aborted.'])
        end
    end
    pool.delete()

    %% Solution
    theta_reg_sol = parameters_reg_sol(:,1:end-5);
    H_isotropic_reg_sol = parameters_reg_sol(:,end-4:end-2);
    H_kinematic_reg_sol = parameters_reg_sol(:,end-1:end);
    theta_save = [theta_save; theta_reg_sol];
    H_isotropic_save = [H_isotropic_save; H_isotropic_reg_sol];
    H_kinematic_save = [H_kinematic_save; H_kinematic_reg_sol];
    cost_save = [cost_save; cost_reg_sol];

else
    error('This method is not implemented.')
end
toc

%% Save Solution
solution.theta = theta_save;
solution.H_isotropic = H_isotropic_save;
solution.H_kinematic = H_kinematic_save;
solution.cost = cost_save;
save(in.save_path,'solution')

end