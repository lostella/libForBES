close all;
clear;

rng(0, 'twister');

% customize the path in the next line
n = 5000;
m = 20000;

w = sprandn(n, 1, 0.3);  % N(0,1), 30% sparse
v = randn(1);            % random intercept

X = sprandn(m, n, 10/n);
btrue = sign(X*w + v);

% noise is function of problem size use 0.1 for large problem
b = sign(X*w + v + sqrt(0.1)*randn(m,1)); % labels with noise

A = spdiags(b, 0, m, m) * [X, ones(m, 1)];

ratio = sum(b == 1)/(m);
lam = 0.1 * norm((1-ratio)*sum(A(b==1,:),1) + ratio*sum(A(b==-1,:),1), 'inf');


fprintf('%d instances, %d features, nnz(A) = %d\n', size(A, 1), size(A, 2), nnz(A));

f = quadLoss(lam);
g = hingeLoss(1, b);
constr = {A, -1, zeros(m,1)};
opt.tol = 1e-10;

opt.method = 'cg-dyhs';
tic; out1 = forbes(f, g, zeros(m, 1), [], constr, opt); toc
out1

opt.method = 'lbfgs';
tic; out2 = forbes(f, g, zeros(m, 1), [], constr, opt); toc
out2

opt.method = 'fbs';
tic; out3 = forbes(f, g, zeros(m, 1), [], constr, opt); toc
out3
