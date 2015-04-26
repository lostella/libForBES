close all;
clear;

% rng(0, 'twister'); % uncomment this to control the random number generator

m = 3000; % number of observations
n = 50000; % number of features
x_orig = sprandn(n, 1, 100/n); % generate random sparse model
A = sprandn(m, n, 500/n); % generate random sparse design matrix
b = A*x_orig + randn(m, 1)/10; % compute labels and add noise

fprintf('%d nonzero features\n', nnz(A));
fprintf('%.2f nnz per row\n', nnz(A)/numel(A)*n);

% for lam >= lam_max the solution is zero
lam_max = norm(A'*b,'inf');
lam = 0.1*lam_max;

% setup problem
prob.A = A;
prob.b = b;
prob.Q = 1;
prob.g = @(x, gam) L1Norm(x, gam, lam);
prob.x0 = zeros(n, 1);

opt.tolOpt = 1e-6;
opt.method = 'lbfgs';
tic; out = minfbe(prob, opt); toc
out
opt.method = 'cg-dyhs';
tic; out = minfbe(prob, opt); toc
out

