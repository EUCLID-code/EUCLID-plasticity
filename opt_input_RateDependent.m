function in = opt_input_RateDependent()
% opt_input_RateDependent provides elementary definitions needed for the optimization problem
%
% ## Comments
%
% Please make sure that the correct path to the FEM data is specified.
%
% ## Input Arguments
%
% _none_
%
% ## Output Arguments
%
% `in` (_struct_) - contains information needed for the optimization
% problem
%
% `in.use_mex` (_logical_) - specifies if .mex file is used for the cost
% evaluation
%
% `in.n_workers` (_double_) - number of workers used for the parallel
% optimization
%
% `in.datafile` (_char_) - filename of the FEM dataset
%
% `in.noise_level` (_char_) - noise level to be considered (e.g., `in.noise_level=''`, `in.noise_level='1'`,
% `in.noise_level='3'`, `in.noise_level='5'`)
%
% `in.filter` (_char_) - temporal data smoothing method (e.g., `in.filter=''` for no
% smoothing, `in.filter='sgolay'` for Savitzky-Golay smoothing)
%
% `in.filter_length` (_double_) - specifies the moving window length for
% the temporal smoothing
%
% `in.library` (_obj_) - model library
%
% `in.q_case` (_double_) - specifies which norm to use during cost function
% calculation (e.g. `in.q_case=0` for L2-norm)
%
% `in.q` (_double_) - specifies which norm to use during cost function
% calculation (e.g. `in.q=2` for L2-norm)
%
% `in.lambda_r` (_double_) - weighting factor for net reaction force
% contributions
%
% `in.lambda_p` (_double_) - weighting factor for the sparsity promoting
% regularization term
%
% `in.p` (_double_) - specifies which norm to use for the sparsity
% promoting regularization
%
% `in.lambda_p_factor` (_double_) - lambda_p is multiplied by this factor
% to obtain solutions for different lambda_p
%
% `in.n_lambda_p` (_double_) - number of different values of lambda_p to be
% considered
%
% `in.n_eval` (_double_) - evaluate cost for every n_eval-th load step
%
% `in.method` (_char_) - method used for optimization
%
% `in.preconditioner` (_char_) - preconditioner used for optimization
% 
% `in.preconditioner_2` (_char_) - preconditioner used for optimization
%
% `in.batch_ratio` (_char_) - ratio defining how to split the data into two batches
%
% `in.MaxFunEvals` (_double_) - maximum minimization function evaluations
% during each optimization process
%
% `in.TolFun` (_double_) - stopping criterion for optimization process
% based on minimization function value
%
% `in.TolX` (_double_) - stopping criterion for optimization process based
% on solution values
%
% `in.MaxTime` (_double_) - maximimum runtime of an optimization process
%
% `in.theta_init` (_double_) - initial guess for the parameters
%
% `in.n_guess` (_double_) - number of random initial guesses
%
% `in.perturbation_guess_theta` (_double_) - specifies the magnitude of
% perturbation applied to the preconditioning solution for generating
% random initial guesses
%
% `in.save_path` (_char_) - path for saving results
%

in.use_mex = false;
in.n_workers = 1;

%% Data
pathname = 'C:\Users\mflas\Documents\Development\DATA_PUBLIC\FEM_data_GSM'; % path to the FEM data
addpath(genpath(pathname));
dataID = 4;
switch dataID
    case 1
        in.datafile = '20220907T130723_DATA_lib_VEVPHardMixLin_LE_plate_elliptic_holes';
    case 2
        in.datafile = '20220907T132348_DATA_lib_VEVPHardMixLin_VE_plate_elliptic_holes';
    case 3
        in.datafile = '20220907T135253_DATA_lib_VEVPHardMixLin_VEEPHardIsoLin_plate_elliptic_holes';
    case 4
        in.datafile = '20220907T134236_DATA_lib_VEVPHardMixLin_LEVPHardKinLin_plate_elliptic_holes';
    case 5
        in.datafile = '20220907T134704_DATA_lib_VEVPHardMixLin_VEVPHardMixLin_plate_elliptic_holes';
end

%% Noise
% in.noise_level = '';
in.noise_level = '1';
% in.noise_level = '3';
% in.noise_level = '5';

if ~isempty(in.noise_level) && dataID > 0
    in.datafile = [in.datafile '_noise_' in.noise_level];
end

%% Data Smoothing
in.filter = '';
in.filter = 'sgolay';
in.filter_length = 10;

if ~isempty(in.noise_level) && ~isempty(in.filter) && dataID > 0
    in.datafile = [in.datafile '_' in.filter '_' num2str(in.filter_length)];
end

%% (Hyper-)Parameters
%% Library
n_Maxwell = [1 1];
in.library = lib_VEVPHardMixLin_v2(n_Maxwell);

%% Overdetermined System of Equations
% case 0: L2-norm (Lq-norm)
% case 1 & 2: Inf-norm
in.q_case = 0;
in.q = 2;
in.lambda_r = 100;

%% Regularization
in.lambda_p = 1e-4;
in.p = 1;
in.lambda_p_factor = 2;
in.n_lambda_p = in.n_workers;

%% Miscellaneous
in.n_eval = 1; % evaluate cost for every n-th load step

%% Optimization Method
in.method = 'lsqnonlin_VEVP_par_perturbation';

%% Preconditioner
in.preconditioner = 'linear_elastic_all_steps';
in.preconditioner_2 = 'none';
in.batch_ratio = 1/2;

%% Miscellaneous
in.MaxFunEvals = 1000;
in.TolFun = 1e-6;
in.TolX = 1e-6;
in.MaxTime = 3; % hours

%% Initial Guesses
in.theta_init = in.library.theta_min;
in.theta_init(1) = 1; % must not be zero
in.theta_init(2) = 1; % must not be zero
in.theta_init(end) = 100; % must not be zero
in.n_guess = in.n_workers;
in.perturbation_guess_theta = ones(size(in.theta_init));

%% Save Path
if contains(in.datafile,'test')
    in.save_path = in.datafile;
else
    in.save_path = in.datafile(22:end);
end
in.save_path = [datestr(now,'yyyymmTddHHMMSS') '_RESULTS_' in.method '_' in.save_path];
if not(isfolder('EUCLID_results'))
    mkdir('EUCLID_results')
end
in.save_path = ['EUCLID_results/' in.save_path];

end