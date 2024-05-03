function[nll_t] = cal_nll_proxy_j(bj0, levels_pi, j, data_order_pi, deltaj, xj)


    n = size(data_order_pi, 1);
    if nargin == 5
        [neg_log_like] = cal_nll_j(bj0, levels_pi, j, data_order_pi);
    else
        [neg_log_like] = cal_nll_j_ia(bj0, levels_pi, j, data_order_pi, xj);
    end
    nll_t = neg_log_like /n  -  dot(cell2mat(deltaj), cell2mat(bj0));
end
