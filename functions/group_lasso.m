function penalty = group_lasso(bj_coefs, levels_pi, j)
dj = length(bj_coefs);

start = 1;  % to separate parents (one hot coding)
penalty = 0;
for i = 1:(j-1)
    di = levels_pi(i) - 1;
    bji = zeros(dj, di);
    % find bji
    for l = 1:dj
        bji(l,:) = bj_coefs{l}(start + 1 : start + di);
    end
    start = start + di;

    bji_norm = norm(reshape(bji, 1, []), 2);
    penalty = penalty + bji_norm;
end

end