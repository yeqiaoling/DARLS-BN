function beta_cur = swap_beta(beta_cur, swap_nodes)

start = 1;
last = length(swap_nodes);
while start < last;
    node1 = swap_nodes(start);
    node2 = swap_nodes(last);
    [beta_cur{node1}, beta_cur{node2}] = deal(beta_cur{node2}, beta_cur{node1});
    start = start + 1;
    last = last - 1;
end

end