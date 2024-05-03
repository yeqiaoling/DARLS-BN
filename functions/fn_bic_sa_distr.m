function [adj_hat_l2_norm_org] = fn_bic_sa_distr(X, varargin)
%
% update fn_bic_sa_distr with choosing internal parameters autamatically
%
% SA: each step of SimulatedAnnealinng is exact solution, which samples a permutation to
%     minimize objective function; 
% SA inputs: 
%
%     min_prop, lambda_num: in BIC path, inidicating solution path (min
%                             lambda and # of lambda)
%    
%% inputs
parser = inputParser;
parser.KeepUnmatched = true;
addOptional(parser, 'version', 3);
addOptional(parser, 'fix_lambda', 0);
% distributed 
addOptional(parser, 'workers', 10);
% bic
addOptional(parser, 'min_prop', 1e-1) 
addOptional(parser, 'lambda_num', 10) 
addOptional(parser, 'BIC', 1)
addOptional(parser, 'use_cd_lambda', 0)
addOptional(parser, 'gamma_rec', -1) 
% SA parameters
addOptional(parser, 'HIGHTEMP', 0);
addOptional(parser,'T_min', 0.00005);
addOptional(parser,'T_max', 0.05);
addOptional(parser,'N', 1e4);
addOptional(parser, 'k', 10);
addOptional(parser, 'step', 0.999);
% getL
addOptional(parser, 'shrinkage', 0.8) 
addOptional(parser, 'TMAX', 50) 
addOptional(parser, 'TOL', 0.05) 

parse(parser, varargin{:});
% distributed 
workers = parser.Results.workers;
% bic
min_prop = parser.Results.min_prop;
lambda_num = parser.Results.lambda_num;
gamma_rec = parser.Results.gamma_rec;
% SA parameters
N =  parser.Results.N;
HIGHTEMP =  parser.Results.HIGHTEMP;
T_max = parser.Results.T_max;
T_min = parser.Results.T_min;
k =  parser.Results.k;
step = parser.Results.step;
% getL
shrinkage =  parser.Results.shrinkage;
TMAX =  parser.Results.TMAX;
TOL =  parser.Results.TOL;

%% initialization
% data 
var_names = X.Properties.VariableNames;
X = table2array(X); [n,p] = size(X);
[X, levels] = factor_to_numeric(X);

mach_ind = randsample(workers, n, true);

% initialize 

ini_topo_sort = randperm(p);

fprintf('Start BIC selection\n')

[~, gamma, lambda] = bic_sel_logistic(X, mach_ind, levels, ini_topo_sort, ...
    'min_prop', min_prop, 'lambda_num', lambda_num,...
    'shrinkage', shrinkage, 'TMAX', TMAX, 'gamma_rec', gamma_rec, ...
    'TOL', TOL);

fprintf('gamma: %2.2f, lambda: %2.4f\n', gamma, lambda)

fprintf('Start SA\n')
rand_seed = randi(1e5);
rng(rand_seed);

% SA
[fvals, ~, ~, pi_hat, beta_hat, ~, ~] = ...
    sa_update(X, levels, mach_ind, ini_topo_sort, ...
    gamma, lambda, 'T_max', T_max, 'T_min',  T_min, 'N', N, 'k', k, ...
    'step', step, 'HIGHTEMP', HIGHTEMP, ...
    'shrinkage', shrinkage, 'TMAX', TMAX, 'TOL', TOL, 'runtime', 0);
[adj_hat_cell, adj_hat_l2_norm, adj_hat_inf_norm, adj_hat] = beta_to_bji(beta_hat, levels(pi_hat));

[adj_hat_org] = convert_to_org_label(pi_hat, adj_hat);
[adj_hat_l2_norm_org] = convert_to_org_label(pi_hat, adj_hat_l2_norm);
[adj_hat_inf_norm_org] = convert_to_org_label(pi_hat, adj_hat_inf_norm);

end
