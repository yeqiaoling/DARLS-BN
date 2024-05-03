function[X_new, levels] = factor_to_numeric(X)
% find the levels of (X_1, ..., X_p), and conver the factor matrix to
%   numeric matrix if applicable
% input: data matrix X
% output: levels of (X_1, ..., X_p)

[n, p] = size(X);
levels = zeros(p, 1);
if ~isnumeric(X(1,1))
    X_new = zeros(n, p);
else 
    X_new = X;
end

for j = 1:p
    vector = X(:,j);
    levels(j) = length(unique(vector));
    if ~isnumeric(vector(1))
        X_new(:,j) = grp2idx(vector);
    end
end

end