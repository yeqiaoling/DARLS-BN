function [] = write_file(array, filename, var_names)
    switch nargin
        case 3
            tab = array2table(array, 'VariableNames', var_names);
            writetable(tab, filename,  'Delimiter', ' ') 
        case 2
            tab = array2table(array);
            writetable(tab, filename,  'Delimiter', ' ', 'WriteVariableNames', false) 
    end
end