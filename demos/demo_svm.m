close all;
clear;

% rng(0, 'twister');

% customize the path in the next line
n = 5000;
m = 20000;

w = sprandn(n, 1, 0.1);  % N(0,1), 10% sparse
v = randn(1);            % random intercept

X = sprandn(m, n, 10/n);
btrue = sign(X*w + v);

% noise is function of problem size use 0.1 for large problem
b = sign(X*w + v + sqrt(0.1)*randn(m,1)); % labels with noise

A = spdiags(b, 0, m, m) * X;

ratio = sum(b == 1)/(m);
lam = 0.1 * norm((1-ratio)*sum(A(b==1,:),1) + ratio*sum(A(b==-1,:),1), 'inf');


fprintf('%d instances, %d features, nnz(A) = %d\n', size(A, 1), size(A, 2), nnz(A));

prob.f = squaredNorm(lam);
prob.g = hingeLoss(1, b);
prob.A = A;
prob.B = -1;
prob.b = zeros(m,1);
% run forbes
opt.tol = 1e-12;
tic;out = forbes(prob, opt);toc
% run acceleated dual proximal gradient
opt_amm.method = 'fbs';
opt_amm.tol = 1e-12;
tic;out_amm = forbes(prob,opt_amm);toc
