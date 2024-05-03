function [beta0] = initialize_beta(levels, varargin)
%% parameters
parser = inputParser;
parser.KeepUnmatched = true;
addOptional(parser,'ia_nums', 10e10) 

parse(parser, varargin{:});
ia_nums =  parser.Results.ia_nums;

% initialize beta: cell(p,1), beta{j}{l}: (1 + d1 + ... + dp), l = 1~dj
% if with interaction terms:  
%       beta{j}{l}: (1 + 2 * (d1 + ... + d_{j-1}) + |ia]
p = length(levels);
beta0 = cell(1, p);

if sum(reshape(ia_nums, 1, [])) == 10e10 && length(ia_nums) == 1
    % each cell: d_j * d_i
    dim2 = sum(levels) - p + 1;
    for j = 1:p
        dj = levels(j) - 1;
        beta0{j} = cell(1, dj);
        for l = 1:dj
            beta0{j}{l} = zeros(1, dim2);
        end
    end
else
    ds = levels - 1;
    for j = 1:p
        dj = levels(j) - 1;
        beta0{j} = cell(1, dj);
        dim2 = 1 + sum(ds(1:max(j-1, 0))) + ia_nums(j);
        for l = 1:dj
            beta0{j}{l} = zeros(1, dim2);
        end
    end
    
end