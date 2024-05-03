function [bic, gamma, lambda] = bic_sel_logistic(X, mach_ind, levels, ini_topo_sort, varargin)
%
% objective: given a permutation, use BIC to select a lambda & gamma 
%           in MCP penalty. At each (gamma, lambda), use regularized 
%           regression to estimate the loss by the proximal gradient 
%           algorithm. 
%           lambda: max --> zero coef; min --> 0,1 * max 
%           gamma: {2, 10, 50, 100}
% inputs: X - design matrix
%         mach_ind - data location indicator
%         var_vat - determine categorical or numerical variable
%         topo_sort - an initial order (topological sort)
% output: bic - minimium bic value
%         (gamma, lambda) - a pair of parameter with smallest BIC, 
%
% BIC = c ln(max(n,p))*k + 2 * neg-log-likelihood
%       where c is a constant with default value of 
%             k is the # of estimates (including variance estimate)
%
%% parameters
parser = inputParser;
parser.KeepUnmatched = true;
addOptional(parser, 'min_prop', 1e-1) 
addOptional(parser, 'lambda_num', 10) 
addOptional(parser, 'c', 1) 
addOptional(parser, 'gamma_rec', -1) 
% getL fn
addOptional(parser, 'shrinkage', 0.8) 
addOptional(parser, 'TMAX', 50) 
addOptional(parser, 'TOL', 1e-1)

parse(parser, varargin{:});
min_prop =  parser.Results.min_prop;
lambda_num =  parser.Results.lambda_num;
gamma_rec = parser.Results.gamma_rec;
% getL fn
shrinkage = parser.Results.shrinkage;
TMAX = parser.Results.TMAX;
TOL = parser.Results.TOL;
%% 
% handle functions
% truncat_coef = @(Bsa)  Bsa .* (abs(Bsa)>  0.1) ;

[n, ~] = size(X);
p = length(levels);
% lambda_max = 1/sqrt(n); gaussian
% lambda_max = find_max_lambda(X, levels, mach_ind, ini_topo_sort);
lambda_max = 0.1;
lambda_min = lambda_max * min_prop;
% lambda_min = 0.01;
lambda_log_rec = log(lambda_min):(log(lambda_max) - log(lambda_min)) ./ (lambda_num-1) ...
    : log(lambda_max);
lambda_rec = exp(lambda_log_rec);
% lambda_rec = [0.01, 0.012, 0.015, 0.02, 0.05, 0.06];
% lambda_rec = lambda_min:(lambda_max - lambda_min)/(lambda_num-1):lambda_max;
gamma_num = length(gamma_rec);
lambda_num = length(lambda_rec);
fprintf('maximum lambda: %2.2f\n', max(lambda_rec));

[beta0] = initialize_beta(levels(ini_topo_sort));
% bicmin_prop 
bic_rec = zeros(gamma_num, lambda_num);
lvals = zeros(gamma_num, lambda_num);
penalties = zeros(gamma_num, lambda_num);
for j = 1:gamma_num  
    gamma = gamma_rec(j);
    for i = 1:lambda_num
        lambda = lambda_rec(i);
        [lval, ~, beta1] = d_getLpg(X, mach_ind, levels, ini_topo_sort, ...
            gamma, lambda, beta0, 1:p, ...
            'shrinkage', shrinkage, 'TMAX', TMAX, 'TOL', TOL);
        lvals(j, i) = lval;
        penalties(j, i) = dim_ml(beta1, levels(ini_topo_sort));
        bic_rec(j, i) = log(max(n,p)) * dim_ml(beta1, levels(ini_topo_sort)) ...
            + 2 * n * lval;
        % bic_rec(j, i) = c * log(max(n,p)) * (nnz(beta1)-p+1) + 2 * n * lval;
        % bic_rec(j, i) = c * log(max(n,p)) * (nnz(beta1)-p+1) + 2 * lval;
    end
end

[bic, ~] = min(bic_rec(:));
[gamma_ind, lambda_ind] = find(bic_rec == bic);
gamma = gamma_rec(gamma_ind(1));
lambda = lambda_rec(lambda_ind(1));

