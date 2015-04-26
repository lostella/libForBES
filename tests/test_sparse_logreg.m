close all;
clear;

% rng(0, 'twister');

m = 6000;
n = 50000;
x_orig = sprandn(n, 1, 100/n);
A = sprandn(m, n, 200/n);
b = 2*(rand(m,1) <= 1./(1+exp(-A*x_orig))) - 1;
[m, n] = size(A);
lam_max = norm(0.5*(A'*b),'inf')/m;
lam = 0.5*lam_max;

prob.C = diag(sparse(b))*A;
prob.f2 = @(w) LogReg(w);
prob.Lf2 = 1/m;
prob.g = @(x, gam) L1Norm(x, gam, lam);
prob.x0 = zeros(n, 1);

opt.tolOpt = 1e-6;
opt.method = 'lbfgs';
tic; out = minfbe(prob, opt); toc
out
opt.method = 'cg-dyhs';
tic; out = minfbe(prob, opt); toc
out

