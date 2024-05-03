function bn = beta_n(b, n, levels)
% calculate b * n, where b is local estimate and n is local sample size
p = length(b);
bn = initialize_beta(levels);

for j = 1:p
    dj = levels(j) - 1;
    for l = 1:dj
        bn{j}{l} = b{j}{l} * n;
    end
end
end