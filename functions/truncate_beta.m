function [bhat_tr] = truncate_beta(bhat, varargin)

% inputs
parser = inputParser;
parser.KeepUnmatched = true;
addOptional(parser,'truncate_prop', 1e-2)

parse(parser, varargin{:});
truncate_prop = parser.Results.truncate_prop;

% truncate bhat
bhat_tr = bhat;
truncate_val = truncate_prop * max(max(abs(bhat_tr)));
bhat_tr(abs(bhat_tr) < truncate_val) = 0;

end