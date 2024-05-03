addpath([pwd '/functions']);
cd(['/data'])

X = readtable('asia_binary.txt'); 
worker = 5;

bhat = fn_bic_sa_distr(X, 'workers', worker, 'N', 100);
bhat_tr = truncate_beta(bhat, 'truncate_prop', 0.1);
bhat_tr