function [bji] = prox_group_lasso_bji(bji, t)
    [dj, di] = size(bji);
    bji_norm = norm(reshape(bji, 1, []), 2);

    bji = max(0, 1 - t / bji_norm) * bji; 
end