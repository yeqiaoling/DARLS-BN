function [] = fn_truncate(bhat, dir_ini, varargin)
% truncate the estimate by threshold 
% inputs
parser = inputParser;
parser.KeepUnmatched = true;
addOptional(parser,'truncate_prop', 1e-2)

parse(parser, varargin{:});
truncate_prop = parser.Results.truncate_prop;

%% find estimates directory

var_names =  bhat.Properties.VariableNames;
bhat = table2array(bhat);

for prop = truncate_prop
    [bhat_tr] = truncate_beta(bhat, 'truncate_prop', prop);
    write_file(bhat_tr, sprintf('%s/beta_tr_%1.2f.txt', dir_ini, prop), var_names)
end


end