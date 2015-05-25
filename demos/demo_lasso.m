close all;
clear;

% rng(0, 'twister'); % uncomment this to control the random number generator

m = 300; % number of observations
n = 5000; % number of features
x_orig = sprandn(n, 1, 20/n); % generate random sparse model
A = sprandn(m, n, 40/n); % generate random sparse design matrix
b = A*x_orig + randn(m, 1)/10; % compute labels and add noise

fprintf('%d nonzero features\n', nnz(A));
fprintf('%.2f nnz per row\n', nnz(A)/numel(A)*n);

% for lam >= lam_max the solution is zero
lam_max = norm(A'*b,'inf');
lam = 0.1*lam_max;

% setup problem
prob.C = A;
prob.d = b;
prob.f = quadratic();
prob.g = l1Norm(lam);
prob.x0 = zeros(n, 1);

opt.display = 1;
opt.maxit = 1000;
opt.tolOpt = 1e-8;
opt.method = 'lbfgs';
tic; out = forbes(prob, opt); toc
out
opt.method = 'cg-dyhs';
tic; out = forbes(prob, opt); toc
out
