function [adj_org] = convert_to_org_label(pi_hat, adj_hat)
    adj_org = adj_hat;
    [~, indx] = sort(pi_hat);
    adj_org = adj_org(indx, indx);
end