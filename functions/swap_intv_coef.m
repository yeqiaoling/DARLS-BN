function [betak] = swap_intv_coef(betak, levels, topo_sort, swap_intv)

if length(swap_intv) == length(levels)
    return
end

% find dk
dk = length(betak);
% find previous permutation
topo_sort_cur = topo_sort;  
topo_sort_cur(swap_intv) = flip(topo_sort_cur(swap_intv));
levels_pi_cur = levels(topo_sort_cur);

left = sum(levels_pi_cur(1: swap_intv(1) - 1) - 1) + 1 + 1;
right = sum(levels_pi_cur(1: swap_intv(length(swap_intv))) - 1) + 1;
for l = 1:dk
    flip_coef = [];
    for j = swap_intv
        start_ind = sum(levels_pi_cur(1:j - 1) - 1) + 1 ;
        dj = levels_pi_cur(j) - 1;
        flip_coef = [betak{l}(start_ind + 1 : start_ind + dj) flip_coef];
    end
    betak{l}(left: right) = flip_coef;
end