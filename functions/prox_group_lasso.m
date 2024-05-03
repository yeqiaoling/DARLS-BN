function [new_z] = prox_group_lasso(zj_lasso_cell, lambda1, t, levels_pi, j)
t = t * lambda1;
new_z = zj_lasso_cell;
dj = levels_pi(j) - 1;

start = 1;
for i = 1:(j-1)
    di = levels_pi(i) - 1;
    bji = zeros(dj, di);
    % find bji: dj * di
    for l = 1:dj
        bji(l,:) = zj_lasso_cell{l}(start + 1 : start + di);
    end
    % applt prox group lasso
    bji = prox_group_lasso_bji(bji, t);
    % convert to the original
    for l = 1:dj
        new_z{l}(start + 1 : start + di) = bji(l, :);
    end
    start = start + di;
end


end
