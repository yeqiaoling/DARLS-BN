function norms = norm_bj(b_diff)

p = length(b_diff);
norms = zeros(1, p);
for j = 1:p
    norms(1, j) = norm(cell2mat(b_diff{j}), 2);
end



end