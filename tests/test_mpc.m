close all;
clear;

% Uncomment the following to control the random number generation
% rng(0);

N = 50; % prediction horizon
n = 12; % no. states
m = 4; % no. inputs

%% ---------------------------------------------------- GENERATE PROBLEM

% Dynamics
A = randn(n, n);
A = 1.15*(A/norm(A));
B = randn(n, m);

% Cost matrices
cond_Q = 1e-2;
cond_R = 1e-2;
Q = diag([rand(n-1,1)+cond_Q; cond_Q]);
R = diag([rand(m-1,1)+cond_R; cond_R]);
[~, Q_f] = dlqr(A,B,Q,R);

% Constraints (boxes)
F = [-speye(n); speye(n); sparse(2*m, n)];
G = [sparse(2*n, m); -speye(m); speye(m)];
xmin = -ones(n, 1); xmax = ones(n, 1);
umin = -0.2*ones(m, 1); umax = 0.2*ones(m, 1);
xminf = -ones(n, 1); xmaxf = ones(n, 1);
c = [-xmin; xmax;  -umin; umax];
F_f = [-speye(n); speye(n)];
c_f = [-xminf; xmaxf];

weights = +inf; % > 0 for soft constraints, +inf for hard constraints

% Setup Hessian (for other solvers to be used)
HessBlocks = cell(2*N+1,1);
for i=1:N
    HessBlocks{2*i-1} = sparse(Q);
    HessBlocks{2*i} = sparse(R);
end
HessBlocks{2*N+1} = Q_f;
Hess = blkdiag(HessBlocks{:});

% Make big (sparse) constraint matrix blocks
FG = sparse([F, G]);
constrBlocks = cell(N+1,1);
for i=1:N
    constrBlocks{i} = FG;
end
constrBlocks{N+1} = sparse(F_f);
C = blkdiag(constrBlocks{:});
d = [repmat(c, N, 1); c_f];

% Setup dynamics constraints (except the initial one)
Aaff = sparse(n+N*n,n+(N-1)*(n+m));
baff = zeros(n+N*n,1);
for i=1:N
    Aaff(i*n+1:(i+1)*n, (i-1)*(n+m)+1:n+i*(n+m)) = [-A, -B, speye(n)];
end

% Initial point (solve LP to choose that)
v = randn(N*(n+m)+n,1);
[some_xu, fval_x0, exitflag_x0] = linprog(v, [], [], Aaff(n+1:end,:), baff(n+1:end,1), [repmat([xmin; umin], N, 1); xminf]+(1e-2), [repmat([xmax; umax], N, 1); xmaxf]-(1e-2));
x0 = some_xu(1:n,1);

% Initial constraint
Aaff(1:n,1:n) = speye(n);
baff(1:n,1) = x0; 

%% ---------------------------------------------------- QUADPROG (interior point)

% Solving with hard constraints only
t0 = tic;
[xu_quadprog, fval, exitflag, out_quadprog] = ...
    quadprog(Hess, zeros(size(Hess,2),1), [], [], Aaff, baff, ...
    [repmat([xmin; umin], N, 1); xminf], [repmat([xmax; umax], N, 1); xmaxf], []);
time_quadprog = toc(t0);
fprintf('\nQUADPROG\n');
fprintf('%15s: %7.4e seconds\n', 'time', time_quadprog);
fprintf('%15s: %7.4e\n', 'obj value', fval);
fprintf('%15s: %7.4e\n', 'infeasibility', max(max(0,C*xu_quadprog-d)));

%% ---------------------------------------------------- GUROBI

% if exist('gurobi')
%     model.Q = 0.5*Hess;
%     model.obj = zeros(size(Hess, 1), 1);
%     model.A = Aaff;
%     model.rhs = baff;
%     model.sense = '=';
%     % Solving with hard constraints only
%     model.lb = [repmat([xmin; umin], N, 1); xminf];
%     model.ub = [repmat([xmax; umax], N, 1); xmaxf];
%     t0 = tic;
%     res_gurobi = gurobi(model);
%     time_gurobi = toc(t0);
% end

%% ---------------------------------------------------- ForBES

opt_forbes.display = 0;
opt_forbes.tolOpt = 1e-8; % this means: norm(residual,inf) <= 1e-3
opt_forbes.maxit = 10000;
res_forbes = solve_mpc(Q, R, Q_f, A, B, F, G, F_f, c, c_f, weights, N, x0, opt_forbes);
time_forbes = res_forbes.time;
fprintf('\nForBES\n');
fprintf('%15s: %7.4e seconds\n', 'time', time_forbes);
fprintf('%15s: %d\n', 'iterations', res_forbes.iterations);
fprintf('%15s: %7.4e\n', 'obj value', 0.5*(res_forbes.x'*(Hess*res_forbes.x)));
fprintf('%15s: %7.4e\n', 'infeasibility', max(max(0,C*res_forbes.x-d)));