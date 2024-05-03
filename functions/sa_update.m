function [fvals, lval, fval, pi_hat, beta_hat, iniP_fval, iniP_beta] = ...
    sa_update(X, levels, mach_ind, ini_topo_sort, gamma, lambda, varargin)
%
% object: use the simulated annealing to sample permutation for minimizing  
%           a loss function, which is negative log-likelihood + MCP 
%           penalty (with parameter gamma & lambda).
% detail: Each P propose is flipping an interval of a permutation of length 
%          'k', then calculate the loss for the proposed one and accept it
%          with the Metropolis rate. 
% note: SA update is performed on corresponding columns of coefficient
%       matrix (instead of the whole matrix). Because of that, each
%       update in the SA is a submatrix of a lower triangular matrix.
% penalty: MCP penalty with gamma & lambda (penalty parameters)
% 
% explicit inputs: X - data matrix
%                  ini_topo_sort - initial topo sort
%                  lambda - penalty function parameter
%                  gamma - penalty parameter 
% inexplicit inputs: SA parameters & update parameters (stoping criteria) & save
%                       indicators (SAVE, FLAG etc)
%                   SAVE: 1 ---> save 
% 
%% inputs
parser = inputParser;
parser.KeepUnmatched = true;
% SA
addOptional(parser,'T_min', 0.0001);
addOptional(parser,'T_max', 0.01);
addOptional(parser,'N', 1e4);
addOptional(parser,'step',0.999)
addOptional(parser,'HIGHTEMP',0)
addOptional(parser,'k', 10)
addOptional(parser,'MAX_COUNT', 100)
% getL fn
addOptional(parser, 'shrinkage', 0.8) 
addOptional(parser, 'TMAX', 1e2) 
addOptional(parser, 'TOL', 1e-1)
% TEST
addOptional(parser,'FLAG', 0)
addOptional(parser,'runtime', 0)

parse(parser, varargin{:});
% SA
T_min =  parser.Results.T_min;
T_max =  parser.Results.T_max;
N =  parser.Results.N;
step =  parser.Results.step;
HIGHTEMP =  parser.Results.HIGHTEMP;
k =  parser.Results.k;
MAX_COUNT = parser.Results.MAX_COUNT;
% getL fn
shrinkage = parser.Results.shrinkage;
TMAX = parser.Results.TMAX;
TOL = parser.Results.TOL;
% TEST
FLAG = parser.Results.FLAG;
runtime = parser.Results.runtime;

%% parameters
% rng(1);

[~, p] = size(X); 
% T_max = T_max * p;
if (HIGHTEMP)
    T_max = 20;
end 
% initial value/estimate
beta_initial = initialize_beta(levels(ini_topo_sort));
[ ~, iniP_fval, iniP_beta] = d_getLpg(X, mach_ind, levels, ini_topo_sort, ...
gamma, lambda, beta_initial, 1:p, ...
'shrinkage', shrinkage, 'TMAX', TMAX, 'TOL', TOL);


% initialization
accProb = @(Edel, T) min(1,exp(-Edel/T));

cost_cur = iniP_fval;
beta_cur = iniP_beta;
topo_sort_cur = ini_topo_sort;

% temperature schedule
T = T_max; Tmax = 0;
T_min = min(T_max/10, T_min);
while (Tmax <= 1)
    Tmax = N/(log(T_min) - log(T_max))*log(step); 
    % # of iters per temperature
    if (Tmax > 1)
        break
    end
    step = 0.9*step;
end

% start
iter = 0; fvals = zeros(1, N);
k0 = min(k, p);
weights = ones(1, k0 - 1) * 2;
not_accept_count = 0;

tic;
while T > T_min && not_accept_count < MAX_COUNT
    %fprintf('temperature: %1d\n',T);
    i = 1;
    while i <= Tmax && not_accept_count < MAX_COUNT
        if not_accept_count > 20 && 0
            k0 = max(k0/2, 4);
            weights = ones(1, k0 - 1) * 2;
            not_accept_count = 0;
        end
        k = randsample(2:k0, 1, true, weights);
        start_pos = randsample((p-k+1),1); 
        swap_intv = start_pos:(start_pos+k-1);
        % proposed permutation: t
        % current permutation: c
        topo_sort_prop = topo_sort_cur;  
        topo_sort_prop(swap_intv) = flip(topo_sort_cur(swap_intv));
        beta_new = swap_beta(beta_cur, swap_intv);
        % update 
        [ ~, cost_new, beta_new] = d_getLpg(X, mach_ind, levels, topo_sort_prop, ...
            gamma, lambda, beta_new, swap_intv, ...
            'shrinkage', shrinkage, 'TMAX', TMAX, 'TOL', TOL);
        ap = accProb(cost_new - cost_cur, T);
        ac = (ap > rand(1));
        if(ac)
            topo_sort_cur = topo_sort_prop;
            cost_cur = cost_new;
            beta_cur = beta_new;
            change = 1;
            weights(k-1) = weights(k-1) + 1;
            not_accept_count = 0;
        end
        i = i+1;
        iter = iter + 1;
        fvals(iter) = cost_cur;
        not_accept_count = not_accept_count + 1;
        if FLAG == 1 && ac
            fprintf('Iter: %d, temp: %d, accept: %0.2f, changed: %d, fval: %2.4e\n', iter, T, ap, ac, cost_cur)
            % weights
            % swap_intv
        end
    end
    T = T*step;
end
SA_time = toc;
fvals = nonzeros(fvals);
pi_hat = topo_sort_cur;
beta_hat = beta_cur;
% final aggregate
[ lval, fval, beta_hat] = d_getLpg(X, mach_ind, levels, pi_hat, ...
        gamma, lambda, beta_hat, 1:p, ...
        'shrinkage', shrinkage, 'TMAX', TMAX, 'TOL', TOL);
