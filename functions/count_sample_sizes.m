function local_sample_sizes  = count_sample_sizes(mach_ind)
K = max(mach_ind);
local_sample_sizes = zeros(K, 1); 
for k = 1:K
    local_sample_sizes(k) = sum(mach_ind == k);
end
end