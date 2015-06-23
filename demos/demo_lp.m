close all;
clear;

rng(0,'twister');


m = 1000;
n = 1500;
A = sprandn(m, n, 100/n);
x_star = (rand(n, 1).*(rand(n, 1) < 0.5));
s_star = (rand(n, 1).*(x_star == 0));
y_star = rand(m, 1);
b = A*x_star;
c = A'*y_star+s_star;

% eliminate null rows
normrow = max(abs(A),[],2);
A = A(normrow>0,:);
b = b(normrow>0,1);
[m, n] = size(A);

%% run GUROBI

fprintf('* GUROBI\n');
model.obj = c;
model.A = A;
model.sense = ['='];
model.rhs = b;
params.outputflag = 0; % quiet
% tic;out_gurobi = gurobi(model, params);toc
tic;[x,y,fval, out] = solve_lp(c,A,b);toc
