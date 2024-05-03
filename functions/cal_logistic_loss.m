function[f, l] = cal_logistic_loss(data, levels_per, b_mle, lambda, gamma, vec)

if nargin == 5
    vec = 1:length(levels_per);
end

[n, dim2] = size(data);

l = 0; penalty = 0;

for j = vec
    bj = b_mle{j}; 
    l = l + cal_nll_j(bj, levels_per, j, data) / n;
    [~, bj_coefs] = sep_bj(bj, dim2);
    if gamma < 0 
        penalty = penalty + lambda * group_lasso(bj_coefs, levels_per, j);
    end
    if isnan(l)
        error('Error in calculating the loss.')
    end
end
f = l + penalty;
end