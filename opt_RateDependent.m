function opt_RateDependent()
% opt_RateDependent runs the optimization problem
%
%    ________ ___  ___ _______ ___     ___ ______
%   /  _____//  / /  //  ____//  /    /  //  _   \
%  /  _____//  /_/  //  /___ /  /___ /  //  /_|  |
% /_______//_______//______//______//__//_______/
% Efficient Unsupervised Constitutive Law Identification and Discovery
% https://euclid-code.github.io/
%
% ## Comments
% 
% Please make sure that the correct input is defined in
% `opt_input_RateDependent`. The results are saved in the directory specified in 
% `opt_input_RateDependent`.
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
in = opt_input_RateDependent(); % load input parameters / hyperparameters / settings etc.

%% Load Data
load(in.datafile,'data')
disp(['Material: ' data.material.name])
disp(['Mesh: ' data.mesh.name])
disp(['Boundary Condition: ' data.bc.name])
disp(['Data: ' in.datafile])
theta_min = in.library.theta_min;

%% Optimization Method
options = optimset( ...
    'MaxFunEvals',in.MaxFunEvals, ...
    'TolFun',in.TolFun, ...
    'TolX',in.TolX, ...
    'Display','iter' ... % 'off' or 'iter'
    );

options_par = optimset( ...
    'MaxFunEvals',in.MaxFunEvals, ...
    'TolFun',in.TolFun, ...
    'TolX',in.TolX, ...
    'Display','off' ... % 'off' or 'iter'
    );

t_opt_before = tic;
disp(['Method: ' in.method])
% to solve the optimization problem, different optimization algorithms may
% be used
if strcmp(in.method,'lsqnonlin_VEVP') ...
        || strcmp(in.method,'lsqnonlin_VEVP_par_perturbation')
    %% Objective Function without and with Regularization
    in_without_reg = in;
    in_without_reg.lambda_p = 0;
    % define objective function without regularization term
    fun_vec_without_reg = def_objective_vec_VEVP(in_without_reg,data);

    % define objective function without regularization term considering
    % only the first time step
    data_1st_step = data;
    data_1st_step.results.time = data.results.time(1);
    data_1st_step.results.u = data.results.u(:,1);
    data_1st_step.results.reaction = data.results.reaction(:,1);
    fun_vec_without_reg_1st_step = def_objective_vec_VEVP(in_without_reg,data_1st_step);

    % define objective function without regularization term considering
    % only the first batch of the data
    data_1st_batch = data;
    data_1st_batch.results.time = data.results.time(1:floor(in.batch_ratio*end));
    data_1st_batch.results.u = data.results.u(:,1:floor(in.batch_ratio*end));
    data_1st_batch.results.reaction = data.results.reaction(:,1:floor(in.batch_ratio*end));
    fun_vec_without_reg_1st_batch = def_objective_vec_VEVP(in_without_reg,data_1st_batch);

    % define objective function with regularization term
    fun_vec = def_objective_vec_VEVP(in,data);

    %% Test True Material Parameters
    % note that the true material must be converted into the right format
    % for cost computation
    if strcmp(in.library.name,'lib_VEVPHardMixLin_v1')
        library_true = lib_VEVPHardMixLin_v1([data.material.n_Maxwell_G,data.material.n_Maxwell_K]);
    elseif strcmp(in.library.name,'lib_VEVPHardMixLin_v2')
        library_true = lib_VEVPHardMixLin_v2([data.material.n_Maxwell_G,data.material.n_Maxwell_K]);
    end
    theta_true = library_true.get_theta_extended(in.library.n_Maxwell,data.material.G0,data.material.K0,data.material.Gi,data.material.gi,data.material.Ki,data.material.ki,data.material.H_iso,data.material.H_kin,data.material.eta,data.material.sigma_0);
    disp(' ')
    disp(['True parameters (theta): ' num2str(theta_true)])

    t_cost_before = tic;
    % compute the objective function for the true parameters
    cost_true = fun_vec_without_reg(theta_true);
    cost_true = sum(cost_true.^2);
    t_cost = toc(t_cost_before); % measure time to compute the objective function
    disp(['Cost using the true parameters: ' num2str(cost_true)])
    disp(['Time needed to compute the cost: ' num2str(t_cost)])

    %% Save Solution Array
    solution.in = in;
    solution.theta_true = theta_true;
    solution.cost_true = cost_true;
    if ~contains(in.datafile,'noise')
        solution.in.standard_deviation = 0;
    else
        solution.in.standard_deviation = data.results.standard_deviation;
    end
    
    %% Initial Guess
    theta_init = in.theta_init;
    theta_save = theta_init;
    disp(' ')
    disp(['Initial guess (theta): ' num2str(theta_init)])
    
    t_cost_init_before = tic;
    % compute the objective function for the initial parameter guess
    cost_init = fun_vec_without_reg(theta_init);
    cost_init = sum(cost_init.^2);
    cost_save = cost_init;
    t_cost_init = toc(t_cost_init_before); % measure time to compute the objective function
    disp(['Cost using the initial parameter guess: ' num2str(cost_init)])
    disp(['Time needed to compute the cost: ' num2str(t_cost_init)])

    %% Preconditioner
    disp(' ')
    disp(['Preconditioner: ' in.preconditioner])
    % to accelerate the optimization process, different preconditioning
    % steps may be used
    % the idea is to solve the problem with a reduced set of parameters
    % first and then use the result as the initial condition for the higher
    % dimensional optimization
    t_opt_precondition_before = tic;
    if strcmp(in.preconditioner,'linear_elastic_1st_step')
        % assume an isotropic linear elastic material and calibrate the
        % stiffness given the information at the first time step
        theta_init(1:2) = lsqnonlin(fun_vec_without_reg_1st_step,theta_init(1:2),theta_min(1:2),[],options);
        % compute the objective function for the initial parameter guess
        cost_init = fun_vec_without_reg(theta_init);
        cost_init = sum(cost_init.^2);
    elseif strcmp(in.preconditioner,'linear_elastic_all_steps')
        % assume an isotropic linear elastic material and calibrate the
        % stiffness given the information at all time steps
        theta_init(1:2) = lsqnonlin(fun_vec_without_reg,theta_init(1:2),theta_min(1:2),[],options);
        % compute the objective function for the initial parameter guess
        cost_init = fun_vec_without_reg(theta_init);
        cost_init = sum(cost_init.^2);
    elseif strcmp(in.preconditioner,'linear_viscoelastic_1st_batch')
        % assume an isotropic linear viscoelastic material and calibrate the
        % stiffness given the information of the first batch of the data
        consider = 1 : (2 + 2*sum(in.library.n_Maxwell));
        theta_init(consider) = lsqnonlin(fun_vec_without_reg_1st_batch,theta_init(consider),theta_min(consider),[],options);
        cost_init = fun_vec_without_reg(theta_init);
        cost_init = sum(cost_init.^2);
    end
    theta_save = [theta_save; theta_init];
    cost_save = [cost_save; cost_init];
    t_opt_precondition = toc(t_opt_precondition_before); % measure time to precondition
    disp(['Parameters after preconditioning (theta): ' num2str(theta_save(end,:))])
    disp(['Cost after preconditioning: ' num2str(cost_save(end))])
    disp(['Time needed to precondition: ' num2str(t_opt_precondition)])

    %% Preconditioner with Fixed Parameters
    if strcmp(in.preconditioner_2,'fix_viscoelastic_parameters')
        t_opt_precondition_before = tic;
        consider = 1 : (2 + 2*sum(in.library.n_Maxwell));
        lb = theta_min; lb(consider) = theta_init(consider);
        ub = Inf(size(lb)); ub(consider) = theta_init(consider);
        [theta_init, min_cost] = lsqnonlin(fun_vec_without_reg,theta_init,lb,ub,options);
        theta_save = [theta_save; theta_init];
        cost_save = [cost_save; min_cost];
        t_opt_precondition = toc(t_opt_precondition_before); % measure time to precondition
        disp(['Parameters after preconditioning with fixed viscoelastic parameters (theta): ' num2str(theta_save(end,:))])
        disp(['Cost after preconditioning with fixed viscoelastic parameters: ' num2str(cost_save(end))])
        disp(['Time needed to precondition with fixed viscoelastic parameters: ' num2str(t_opt_precondition)])
    end

    %% Optimization
    if strcmp(in.method,'lsqnonlin_VEVP')
        [theta, min_cost] = lsqnonlin(fun_vec,theta_init,theta_min,[],options);
        theta_save = [theta_save; theta];
        cost_save = [cost_save; min_cost];
    elseif strcmp(in.method,'lsqnonlin_VEVP_par_perturbation')

        if strcmp(in.method,'lsqnonlin_VEVP_par_perturbation')
            %% Random Initial Guesses
            % generate random initial guesses by perturbing the solution obtained
            % from the preconditioning step
            theta_random = repmat(theta_init,in.n_guess,1);
            for idx = 2:in.n_guess
                admissible = false;
                while ~admissible
                    % random perturbations are assumed to be positive to avoid
                    % negative parameters
                    theta_random(idx,:) = theta_init + in.perturbation_guess_theta.*abs(randn(size(theta_init)));
                    % check if random initial guess is admissible
                    if (theta_random(idx,1) > theta_min)
                        admissible = true;
                    end
                end
            end
        elseif strcmp(in.method,'lsqnonlin_VEVP_par_random')
            % Another strategy for generating random initial guesses would be
            % to optimize the unregularized minimization problem to obtain a
            % solution theta_opt and afterwards choose random initial guesses
            % between theta_min and 2*theta_opt. In this way the orders of
            % magnitude of the parameters could be estimated approximately.
        end

        %% Parallel Optimization without Regularization
        % minimize the objective function without regularization term for the
        % randomly generated initial guesses
        % optimizations for different initial guesses are executed in parallel
        theta_random_sol = zeros(size(theta_random));
        cost_random_sol = zeros(in.n_guess,1);
        cluster = parcluster('local');
        pool = parpool(cluster,in.n_workers);

        disp(['Total number of parallel optimizations with different initial guesses:' num2str(in.n_guess)])
        parfor idx = 1:in.n_guess
            try
                [theta, min_cost] = lsqnonlin(fun_vec_without_reg,theta_random(idx,:),theta_min,[],options_par);
                theta_random_sol(idx,:) = theta;
                cost_random_sol(idx) = min_cost;
                disp(['Parallel optimization number ' num2str(idx) ' converged.'])
            catch
                % if lsqnonlin fails to solve the problem, the parameters are
                % set to NaN
                theta_random_sol(idx,:) = NaN;
                cost_random_sol(idx) = Inf;
                disp(['Parallel optimization number ' num2str(idx) ' aborted.'])
            end
        end

        %% Save Solution
        theta_save = [theta_save; theta_random_sol];
        cost_save = [cost_save; cost_random_sol];

        %% Solution Selection
        % select the solution with the lowest objective function value
        % the selected solution will later be used as the initial guess for the
        % regularized optimization
        [~,pos] = min(cost_random_sol);
        theta_init = theta_random_sol(pos,:);

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
        theta_reg_sol = zeros(in.n_lambda_p,length(theta_init));
        cost_reg_sol = zeros(in.n_lambda_p,1);

        disp(['Total number of parallel optimizations with varying regularization:' num2str(in.n_lambda_p)])
        parfor idx = 1:in.n_lambda_p
            in_mod = in;
            in_mod.lambda_p = lambda_p_set(idx);
            fun_vec = def_objective_vec_VEVP(in_mod,data);
            try
                [theta, min_cost] = lsqnonlin(fun_vec,theta_init,theta_min,[],options_par);
                theta_reg_sol(idx,:) = theta;
                cost_reg_sol(idx) = min_cost;
                disp(['Parallel optimization number ' num2str(idx) ' converged.'])
            catch
                theta_reg_sol(idx,:) = NaN;
                cost_reg_sol(idx) = Inf;
                disp(['Parallel optimization number ' num2str(idx) ' aborted.'])
            end
        end
     
        pool.delete()
        
        %% Save Solution
        theta_save = [theta_save; theta_reg_sol];
        cost_save = [cost_save; cost_reg_sol];

    end

    disp('Parameters after optimization (theta): ')
    disp(num2str(theta_save))
    disp('Cost after optimization: ')
    disp(num2str(cost_save))

else
    error('Not implemented')
end
t_opt = toc(t_opt_before);
disp(['Time needed to optimize: ' num2str(t_opt)])

%% Save Solution
solution.theta = theta_save;
solution.cost = cost_save;
save(in.save_path,'solution')

end