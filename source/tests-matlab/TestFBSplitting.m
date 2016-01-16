%% boxqp

Q = [7, 2, -2, -1; 2, 3, 0, -1; -2, 0, 3, -1; -1, -1, -1, 1];
q = [1, 2, 3, 4]';
gam = 0.1;
lb = -1;
ub = +1;
x1 = [+0.5, +1.2, -0.7, -1.1]';
x2 = [-1.0, -1.0, -1.0, -1.0]';
f = quadratic(Q, q);
g = indBox(lb, ub);
opt = [];
opt.method = 'fbs';
opt.variant = 'basic';
opt.maxit = 1000;
opt.tol = 1e-14;
opt.display = 2;
opt.Lf = 1/gam;
opt.adaptive = 0;
out_boxqp = forbes(f, g, x1, [], [], opt);
xstar_boxqp = out_boxqp.x;
% xstar = [ -0.352941176470588,
%           -0.764705882352941,
%           -1.000000000000000,
%           -1.000000000000000]';

%% lasso

A = [1, -2, 3, -4, 5; 2, -1, 0, -1, 3; -1, 0, 4, -3, 2; -1, -1, -1, 1, 3];
b = [1, 2, 3, 4]';
gam = 0.01;
x1 = zeros(5,1);
f = quadLoss(1, b);
g = l1Norm(5);
opt = [];
opt.method = 'fbs';
opt.variant = 'basic';
opt.maxit = 10000;
opt.tol = 1e-14;
opt.display = 2;
opt.Lf = 1/gam;
opt.adaptive = 0;
out_lasso = forbes(f, g, x1, {A}, [], opt);
xstar_lasso = out_lasso.x;
% xstar = [ -0.010238907849511,
%           0,
%           0,
%           0,
%           0.511945392491421]';

%% sparselogreg

A = [1, -2, 3, -4, 5; 2, -1, 0, -1, 3; -1, 0, 4, -3, 2; -1, -1, -1, 1, 3];
b = [1, -1, 1, -1]';
gam = 0.1;
x1 = zeros(5,1);
f = logLoss(1);
g = l1Norm(1);
opt = [];
opt.method = 'fbs';
opt.variant = 'basic';
opt.maxit = 1000;
opt.tol = 1e-14;
opt.display = 2;
opt.Lf = 1/gam;
opt.adaptive = 0;
out_sparselogreg = forbes(f, g, x1, {A, -b}, {}, opt);
xstar_sparselogreg = out_sparselogreg.x;
% xstar = [ 0.0
%           0.0
%           0.215341883018748
%           0.0
%           0.675253988559914]';