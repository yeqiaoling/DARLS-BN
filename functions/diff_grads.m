function [delta] = diff_grads(grad_local, grad, n_local)
% given gradients from local machine, find the weighted average of
% gradients
p = length(grad);
delta = cell(1, p);

for j = 1:p
    dj = length(grad{j});
    delta{j} = cell(1, dj);
    for l = 1:dj
        delta{j}{l} = grad_local{j}{l} ./ n_local - grad{j}{l};  
    end
end
  
end
