function [ l, f, bold] = get_beta_logistic(data, levels, ini_topo_sort, delta, ...
    gamma, lambda, bold, swap_intv, varargin)
%
% objective: given a topological sort, we use proximal gradient algorithm 
% to find a minimizer beta that minimize the MCP regularized logistic
% negative likelihood
%
% inputs: data - design matrix
%         topo_sort - a topolgical sort estimate (proposed)
%         gamma - MCP penalty concavity
%         lambda - MCP penalty strength
%         beta0 - initial beta matrix (p*p) in the order of
%                   topo_sort, suppose propose an interval by 
%                   flipping columns k to m, then only update k to m 
%                   columns in 'beta0'
%         vec - which colomns to updates in 'beta0', i.e. k to m in the
%               previous example
%
% outputs: fval - minimum loss 
%          gval - minimum negative log-likelihood 
%          lres - minimizer
%
%% parameters
parser = inputParser;
parser.KeepUnmatched = true;
addOptional(parser,'shrinkage', 0.8) 
addOptional(parser,'iter_pg', 1e2) 
addOptional(parser,'TOL', 1e-2) 
% testing purpose
addOptional(parser,'FLAG', 0) 

parse(parser, varargin{:});
shrinkage =  parser.Results.shrinkage;
iter_pg =  parser.Results.iter_pg;
TOL =  parser.Results.TOL;
% testing purpose
FLAG =  parser.Results.FLAG;

%%
[n, p] = size(data);
data_order_pi = data(:, ini_topo_sort);
levels_pi = levels(ini_topo_sort);
data_order_pi = one_hot_data(data_order_pi, levels_pi);
data_order_pi = [ones(n, 1) data_order_pi];
[~, dim2] = size(data_order_pi);

% set up
for j = swap_intv
    bj = bold{j}; 
    t = max(1, fro_norm_beta(bj));
    iter = 0;  err = Inf;  
    while iter < iter_pg && err >= TOL; 
        iter = iter + 1;
        bj0 = bj;
        % find gradient
        [tgradj, gradj, inters, valids] = grad_beta_j(bj0, levels_pi, j, data_order_pi, delta{j});
        grad_nm = norm(cell2mat(tgradj),'fro');
        T2 = 0;
        % start line search
        if FLAG
            fprintf('Line search started\n')
        end
        while T2 <= iter_pg
            T2 = T2 + 1;
            bt = cell2mat(bj0) - t / grad_nm * cell2mat(tgradj);
            z = bt;
            % valid z: off-diag
            z = z .* cell2mat(valids);
            % prox_th -- group lasso (NO MCP YET)
            [z_cell] = mat_to_cell(z, levels_pi(j) - 1, dim2);
            [z_cell] = prox_group_lasso(z_cell, lambda, t/grad_nm, levels_pi, j);
            z = cell2mat(z_cell);
            % bt_inters = bt .* (cell2mat(inters) == 1);
            % z_inters = bt_inters;
            % z = cell2mat(z_coef) + z_inters;
            % find loss
            % [z_cell] = mat_to_cell(z, levels_pi(j) - 1, dim2);
            loss_z = cal_nll_proxy_j(z_cell, levels_pi, j, ...
                  data_order_pi, delta{j});
            loss_bold = cal_nll_proxy_j(bj0, levels_pi, j, ...
                data_order_pi, delta{j});
            
            if FLAG
                fprintf('    iteration: %d, step size: %d\n', T2, t)
                fprintf('    LHS: %f, RHS: %f\n', loss_z, ...
                    loss_bold + dot(cell2mat(tgradj), (z - cell2mat(bj0))) ...
                    + 0.5 * norm(z - cell2mat(bj0), 'fro')^2 * grad_nm / t )
            end 
            
            if (loss_z <= loss_bold + dot(cell2mat(tgradj), (z - cell2mat(bj0))) ...
                     + 0.5 * norm(z - cell2mat(bj0), 'fro')^2 * grad_nm / t )
                break
            else
                t = t * shrinkage;
            end
        end
        bj = z_cell;
        err = bj_error(bj, bj0) / max(1, norm(cell2mat(bj),'fro'));
    end
    bold{j} = bj;
end
% calculate loss -- with group lasso penalty
[f, l] = cal_logistic_loss(data_order_pi, levels_pi, bold, lambda, gamma);

