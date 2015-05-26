close all;
clear;

rng(0);

S = system_masses(3);

% Get dynamics from the system
A = S.A;
B = S.B;
[n, m] = size(B);

% Define cost matrices
Q = 1*eye(n);
R = 1*eye(m);
[~, Q_f] = dlqr(A,B,Q,R);

% Input constraints
umin = -10*ones(m,1);
umax = +10*ones(m,1);
xmin = -1*ones(n,1);
xmax = +1*ones(n,1);
% Initial point and prediction horizon
x0 = randn(n, 1);
N = 20;

% Define ForBES problem
prob.f1 = lqrCost(x0, Q, R, Q_f, A, B, N);

% % many blocks
% prob.g = separableSum(  {distBox(xmin,xmax,1e5*ones(n,1)),...
%     indBox(umin,umax), distBox(-ones(n,1),ones(n,1),1e5*ones(n,1))}, ...
%     [repmat([1, 2], 1, N), 3], ...    % indices
%     [repmat([n, m], 1, N), n]);       % dimensions

% % two blocks
% lb = repmat([xmin;umin],N,1);
% ub = repmat([xmax;umax],N,1);
% weights = repmat([1e5*ones(n,1);inf*ones(m,1)],N,1);
% prob.g = separableSum(  {distBox(lb,ub,weights), distBox(-ones(n,1),ones(n,1),1e5*ones(n,1))}, ...
%     [1, 2], ...    % indices
%     [N*(n+m), n]);       % dimensions

% one block
lb = [repmat([xmin;umin],N,1);xmin];
ub = [repmat([xmax;umax],N,1);xmax];
weights = [repmat([1e2*ones(n,1);inf*ones(m,1)],N,1);1e2*ones(n,1)];
prob.g = distBox(lb,ub,weights);

prob.A1 = 1;
prob.B = -1;
prob.b = zeros(n*(N+1)+m*N,1);

% Call solver
opt.display = 1;
opt.tol = 1e-8;
tic;out = forbes(prob,opt);toc
% Call fast AMM
opt_amm.display = 1;
opt_amm.tol = 1e-8;
opt_amm.method = 'fbs';
tic;out_amm = forbes(prob,opt_amm);toc
