function [grad] = average_grad(grads, n)
% given gradients from local machine, find the weighted average of
% gradients
workers = length(grads);
p = length(grads{1});
grad = cell(1, p);

for j = 1:p
    dj = length(grads{1}{j});
    grad{j} = cell(1, dj);
    for l = 1:dj
        cur_sum = 0;
        for worker = 1:workers
            cur_sum = cur_sum + grads{worker}{j}{l};
        end
        grad{j}{l} = cur_sum/n;
    end
end
  
end
