function[err] = bj_error(bj, bold)
dj = length(bj);
err = 0;
for l = 1:dj
    err = max(norm(bj{l} - bold{l}, 2), err); 
end
end
