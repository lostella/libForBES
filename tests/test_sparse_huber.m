close all;
clear;

% rng(0, 'twister');

m = 5000;
n = 20000;
x_orig = sprandn(n, 1, 100/n);
A = sprandn(m, n, 300/n);
outl = rand(m,1) < 0.01;
% add small noise for everyone + large noise for outliers
b = A*x_orig + randn(m, 1)/10 + outl.*randn(m, 1)*10;
lam_max = norm(A'*b,'inf');
lam = 0.1*lam_max;
% since we know what small/large noise means, would do cross-validation otherwise I guess
del = 1;

prob.C = A;
prob.d = b;
prob.f2 = @(x) Huber(x, del);
prob.Lf2 = 1;
prob.g = @(x, gam) L1Norm(x, gam, lam);
prob.x0 = zeros(n, 1);

opt.tolOpt = 1e-6;
opt.method = 'lbfgs';
tic; out = minfbe(prob, opt); toc
out
opt.method = 'cg-dyhs';
tic; out = minfbe(prob, opt); toc
out

