function [frob_norm] = fro_norm_beta(bj) 
% find the frobenius norm of a cell: sqrt(sum of abs(x))
rj = length(bj);
cur_sum = 0;
for l = 1:rj
    cur_sum = cur_sum + norm(bj{l}, 'fro')^2;
end
frob_norm = sqrt(cur_sum);
end
