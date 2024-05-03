function[bj_inters, bj_coefs] = sep_bj(bj, dim2)
% seperate bj to intercepts and coefficients 
%      b/c we penalize coefficients
dj = length(bj);
bj_inters = cell(1, dj);
bj_coefs = cell(1, dj);

for l = 1:dj
    bj_inters{l} = zeros(1, dim2);
    bj_inters{l}(1, 1) = bj{l}(1);
    
    bj_coefs{l} = zeros(1, dim2);
    bj_coefs{l}(1, 2:dim2) = bj{l}(2:dim2);
end
        
end