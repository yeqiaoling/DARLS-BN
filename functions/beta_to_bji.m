function [adj_cell, adj_l2_norm, adj_inf_norm, adj_mat] = beta_to_bji(beta_hat, levels)
% convert beta (cell 1*p) to bji (cell p*p)
% adj_cell: save paramters
% adj_l2_norm: l2 norm
% adj_inf_norm: max norm
% adj_mat: adjacency matrix
% 
p = length(beta_hat);
adj_cell = cell(p, p);
adj_mat = zeros(p, p);
adj_l2_norm = zeros(p, p);
adj_inf_norm = zeros(p, p);
for j = 2:p
    dj = levels(j) - 1;
    start = 1;
    for i = 1:(j-1)
        di = levels(i) - 1;
        bji = zeros(dj, di);
        for l = 1:dj
            bji(l, :) = beta_hat{j}{l}(start + 1: start + di);
        end
        start = start + di;
        adj_cell{j, i} = bji;
        l2_norm = norm(reshape(bji, 1, []),2);
        inf_norm = norm(reshape(bji, 1, []), inf);
        adj_l2_norm(j, i) = l2_norm;
        adj_inf_norm(j, i) = inf_norm;
        adj_mat(j,i) = (l2_norm > 0);
    end
end

end