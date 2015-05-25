close all;
clear;

% rng(0, 'twister'); % uncomment this to control the random number generator

m = 600;
n = 5000;
x_orig = sprandn(n, 1, 30/n);
A = sprandn(m, n, 50/n);
b = 2*(rand(m,1) <= 1./(1+exp(-A*x_orig))) - 1;
lam_max = norm(0.5*(A'*b),'inf')/m;
lam = 0.5*lam_max;

prob.f = logLogistic(1/m);
prob.C = diag(sparse(b))*A;
prob.g = l1Norm(lam);
prob.x0 = zeros(n, 1);

opt.display = 1;
opt.maxit = 1000;
opt.tolOpt = 1e-8;
opt.method = 'lbfgs';
tic; out_lbfgs = forbes(prob, opt); toc
out_lbfgs
opt.method = 'cg-dyhs';
tic; out_cg = forbes(prob, opt); toc
out_cg
opt.method = 'fbs';
tic; out_fbs = forbes(prob, opt); toc
out_fbs
