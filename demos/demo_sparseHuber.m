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

f = huberLoss(del);
g = l1Norm(lam);
opt.display = 1;
opt.tol = 1e-10;

opt.method = 'lbfgs';
out1 = forbes(f, g, zeros(n,1), {A, -b}, [], opt);
out1.ts(end)
out1

opt.method = 'cg-dyhs';
out2 = forbes(f, g, zeros(n,1), {A, -b}, [], opt);
out2.ts(end)
out2

opt.method = 'fbs';
out3 = forbes(f, g, zeros(n,1), {A, -b}, [], opt);
out3.ts(end)
out3
