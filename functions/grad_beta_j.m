function [tgradj, gradj, inters, valids] = grad_beta_j(bj0, levels_pi, j, data_order_pi, deltaj)

[n, ~] = size(data_order_pi);
dj = levels_pi(j) - 1;
gradj = cell(1, dj);
tgradj = cell(1, dj);
inters = cell(1, dj);
valids = cell(1, dj);
muj = zeros(n, dj);

start = 1 + sum(levels_pi(1:(j - 1))) - (j - 1);
xj = data_order_pi(:, start + 1: start + dj); % n * dj
% find mu_jl
for l = 1:dj
    muj(:, l) = exp(data_order_pi * bj0{l}');     
end
row_sum = sum(muj, 2) + 1;
muj = bsxfun(@rdivide, muj, row_sum(:)); % n * dj
for l = 1:dj
    gradj{l} = (xj(:, l) - muj(:,l))' * data_order_pi;
    gradj{l} = - gradj{l} / n;
    tgradj{l} = gradj{l} - deltaj{l};
    
    inters{l} = zeros(size(gradj{l}));
    inters{l}(1, 1) = 1;
    valids{l} = zeros(size(gradj{l}));
    valids{l}(1, 1: start) = 1;
end
    
end