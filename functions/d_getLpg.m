function [ lval, fval, beta1] = d_getLpg(X, mach_ind, levels, ini_topo_sort, ...
    gamma, lambda, beta_initial, swap_intv, varargin)
%
% objective: given a permutation P, use proximal gradient algorithm find 
% a minimizer beta that is compatible with P (beta denotes a DAG with order 
% pi using logistic regression)
%
% loss: negative log-likelihood with MCP penalty 
%
% inputs: X - design matrix 
%         mach_ind - indicates location of data observations
%         topo_sort - a topological sort 9proposed)
%         gamma - MCP penalty concavity
%         lambda - MCP penalty strength
%         beta0 - initial matrix (p+1)*p in the order of
%                 pi
%         vec - which colomns to updates in 'beta0'
%
% outputs: fval - minimum loss 
%          gval - minimum negative log-likelihood 
%         lres - minimizer
%
%% parameters
parser = inputParser;
parser.KeepUnmatched = true;
addOptional(parser,'shrinkage', 0.8) 
addOptional(parser,'TMAX', 1e2) 
addOptional(parser,'TOL', 1e-2) 
addOptional(parser,'iter_pg', 50) 

parse(parser, varargin{:});
shrinkage =  parser.Results.shrinkage;
TMAX =  parser.Results.TMAX;
TOL =  parser.Results.TOL;
iter_pg = parser.Results.iter_pg;

%%
[n, ~] = size(X);
p = length(levels);
workers = max(mach_ind);

iter = 0; err = Inf;  
betas = cell(workers, 1); 
levels_pi = levels(ini_topo_sort);
for worker = 1:workers 
    betas{worker} = initialize_beta(levels_pi);
end
grads = betas;
% betas = zeros(p+1, p, workers);
% grads = zeros(p+1, p, workers);
ns = zeros(workers, 1);
fs = zeros(workers, 1);
ls = zeros(workers, 1);

beta1 = beta_initial;
errs = [];
while iter < TMAX && err >= TOL;
    iter = iter + 1;
    bold = beta1;
    parfor m = 1:workers
        data = X(mach_ind == m, ini_topo_sort);
        [nk, ~] = size(data);
        data = one_hot_data(data, levels_pi);
        data = [ones(nk, 1) data];
        % find gradients
        grads{m} = grad_beta(bold, levels_pi, data);
        ns(m) = nk;
    end
    
    [grad] = average_grad(grads, sum(ns));
    
    parfor m = 1:workers
        data = X(mach_ind == m, :);
        % find delta
        [delta] = diff_grads(grads{m}, grad, ns(m));
        % get beta
        [ ~, ~, b] = get_beta_logistic(data, levels, ini_topo_sort, ...
            delta, gamma, lambda, bold, swap_intv, ...
            'shrinkage', shrinkage, 'iter_pg', iter_pg, 'TOL', TOL);
        % betas{m} = b * ns(m)
        betas{m} = beta_n(b, ns(m), levels_pi);
        % fs(m) = f; ls(m) = l;
    end
    beta1 = average_grad(betas, n);
    % beta1 = sum(betas, 3) / n;
    b_diff = diff_grads(beta1, bold, 1);
    err = max(norm_bj(b_diff) ./ max(1, norm_bj(bold))); 
    % err = max(vecnorm(beta1 - bold)./max(1, norm(bold,'fro')));
    errs(iter) = err;
end

% change indexes of beta_k for k after the swap interval
start = max(swap_intv) + 1;
for k = start:p
    [beta1{k}] = swap_intv_coef(beta1{k}, levels, ini_topo_sort, swap_intv);
end

% calculate loss
parfor m = 1:workers
    data = X(mach_ind == m, ini_topo_sort);
    data = one_hot_data(data, levels_pi);
    data = [ones(ns(m), 1) data];
    [ fs(m), ls(m)] = cal_logistic_loss(data, levels_pi, beta1, lambda, gamma, 1:p);
end
fval = sum(fs .* ns) / n;
lval = sum(ls .* ns) / n;

if 1 == 0
    topo_sort_cur = ini_topo_sort;  
    topo_sort_cur(swap_intv) = flip(topo_sort_cur(swap_intv));
    levels_pi_cur = levels(topo_sort_cur);

    parfor m = 1:workers
        % swap interval
        data = X(mach_ind == m, ini_topo_sort);
        data = one_hot_data(data, levels_pi);
        data = [ones(ns(m), 1) data];
        [f_swap, l_swap] = cal_logistic_loss(data, levels_pi, beta1, lambda, gamma, swap_intv);
        % other elements
        data_pi_cur = X(mach_ind == m, topo_sort_cur);
        data_pi_cur = one_hot_data(data_pi_cur, levels_pi_cur);
        data_pi_cur = [ones(ns(m), 1) data_pi_cur];
        [f_exist, l_exist] = cal_logistic_loss(data_pi_cur, levels_pi_cur, ...
            beta1, lambda, gamma, setdiff(1:length(levels_pi_cur), swap_intv));
        % overall
        fs(m) = f_exist + f_swap;
        ls(m) = l_exist + l_swap;
    end
    fval = sum(fs .* ns) / n;
    lval = sum(ls .* ns) / n;
end
% fprintf('FULL PROX. GRAD. Time: %d, obj: %2.2e, step: %2.2e, error: %2.2e, \n', T, fval, t, err);
% fprintf('g(l): %d, mcp: %2.2e\n', g(l) , MCP(tril(l,-1)));



