function[new_z] = mat_to_cell(z, dj, dim2)

new_z = cell(1, dj);

start = 0;
for l = 1:dj
    new_z{l} = z(start + 1 : start + dim2);
    start = start + dim2;
end
end