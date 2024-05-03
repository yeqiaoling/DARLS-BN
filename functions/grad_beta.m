function [grad] = grad_beta(bold, levels_pi, data) 
% b: a cell (p, p)

p = length(levels_pi);
start = 1;
[n, ~] = size(data);
grad = cell(1, p);
for j = 1:p
    dj = levels_pi(j) - 1;
    grad{j} = cell(1, dj);
    muj = zeros(n, dj);
    xj = data(:, start + 1: start + dj);
    start = start + dj;
    % find mu_jl
    for l = 1:dj
        muj(:, l) = exp(data * bold{j}{l}');     
    end
    row_sum = sum(muj, 2) + 1;
    muj = bsxfun(@rdivide, muj, row_sum(:)); % n * dj
    for l = 1:dj
       grad{j}{l} = -(xj(:,l) - muj(:,l))' * data;
    end
end

end