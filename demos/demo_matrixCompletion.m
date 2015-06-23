% solve a matrix completion problem using forbes

close all;
clear;

rng(0, 'twister'); % uncomment this to control the random number generator

m = 200; % number of rows
n = 200; % number of column of the original matrix M
d = 15000/(m*n); % density of coefficients sampled from M
r = 10; % rank of M

U = randn(m, r);
V = randn(n, r);
M = U*V';

P = sprand(m, n, d) ~= 0; % sampling pattern
B = full(M.*P);

lam = 1e0;

f = quadLoss(1, P(:), B(:));
g = nuclearNorm(m, n, lam);
x0 = zeros(m*n, 1);
opt.tol = 1e-12;

opt.method = 'lbfgs';
tic; out = forbes(f, g, x0, [], [], opt); toc;
out

opt.method = 'fbs';
tic; out = forbes(f, g, x0, [], [], opt); toc;
out
