function[X_dummy] = one_hot_data(X, levels)
% one hot code data X
% input: X_j: r_j level indicated in levels
% output: represent X_j using r_j - 1 dummy variables

[~, p] = size(X);
X_dummy = [];

for j = 1:p 
    dj = levels(j) - 1;
    vector = X(:, j);
    syms = unique(vector);
    for l = 1:dj
        add_col = (vector == syms(l));
        X_dummy = [X_dummy add_col];
    end
end

end