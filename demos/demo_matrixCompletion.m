% solve a matrix completion problem using forbes
close all;
clear;

rng(0, 'twister');

%% ---------------------------------------------------- GENERATE PROBLEM
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

prob.f = quadratic(diag(P(:)));
prob.d = B(:);
prob.g = nuclearNorm(m, n, lam);
prob.x0 = zeros(m*n, 1);
tic; out = forbes(prob); toc
