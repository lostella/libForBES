close all;
clear;

% rng(0, 'twister');

m = 6000;
n = 50000;
x_orig = sprandn(n, 1, 100/n);
A = sprandn(m, n, 200/n);
b = 2*(rand(m,1) <= 1./(1+exp(-A*x_orig))) - 1;
[m, n] = size(A);
lam = 1e0;

prob.x1step = @(w) w/lam;
prob.zstep = @(y, gam) Hinge(y/gam, 1/gam, b);
prob.A1 = A;
prob.B = -1;
prob.c = zeros(m,1);

opt.tolOpt = 1e-6;
opt.method = 'lbfgs';
tic; out = miname(prob, opt); toc
out
opt.method = 'cg-dyhs';
tic; out = miname(prob, opt); toc
out

