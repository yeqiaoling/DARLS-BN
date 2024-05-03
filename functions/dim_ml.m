function dim = dim_ml(b, levels)
dim = 0;
p = length(levels);
for j = 1:p
    nnz(cell2mat(b{j}));
    % dim = dim + nnz(cell2mat(b{j})) - levels(j) + 1;
    dim = dim + nnz(cell2mat(b{j}));
end
end