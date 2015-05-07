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

% pack up problem
[LRs, Ks, Ms, Ls] = RiccatiFactor(Q, R, Q_f, A, B, N);
prob.f1 = lqrCost(x0, A, B, LRs, Ks, Ms, Ls, N);
prob.g = softPositiveOrthant(weights);
prob.A1 = C;
prob.B = 1;
prob.b = [repmat(c, N, 1); c_f];

opt.display = 2;
opt.method = 'lbfgs';
opt.tolOpt = 1e-3;
out = forbes(prob, opt);
opt.method = 'cg-dyhs';
out = forbes(prob, opt);