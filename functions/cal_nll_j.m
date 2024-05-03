function[neg_log_like] = cal_nll_j(bj, levels, j, data)
        
[n, ~] = size(data);
dj = levels(j) - 1;
thetaj = zeros(n, dj);

start = 1 + sum(levels(1:(j - 1))) - (j - 1);
xj = data(:, start + 1 : start + dj);
% find mu_jl

for l = 1:dj
    thetaj(:, l) = data * bj{l}'; % n * 1
end

term1 = sum(sum(thetaj .* xj));
term2 = -sum(log(sum(exp(thetaj), 2) + 1));

neg_log_like = -(term1 + term2);

end