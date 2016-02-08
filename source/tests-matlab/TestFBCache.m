%% boxqp

Q = [7, 2, -2, -1; 2, 3, 0, -1; -2, 0, 3, -1; -1, -1, -1, 1];
q = [1, 2, 3, 4]';
gam1 = 0.1;
gam2 = 0.2;
lb = -1;
ub = +1;
x1 = [+0.5, +1.2, -0.7, -1.1]';
x2 = [-1.0, -1.0, -1.0, -1.0]';
f = quadratic(Q, q);
g = indBox(lb, ub);
proxg = g.makeprox();

gradfx1 = f.Q*x1 + f.q;
fx1 = 0.5*(x1'*(gradfx1 + f.q));
y = x1 - gam1*gradfx1;
[z, gz] = proxg(y, gam1);
diff = x1-z;
FBEx1_gam1 = fx1 + gz - gradfx1'*diff + (0.5/gam1)*(diff'*diff);
gradFBEx1_gam1 = diff/gam1 - f.Q*diff;

y = x1 - gam2*gradfx1;
z = proxg(y, gam2);
diff = x1-z;
FBEx1_gam2 = fx1 + gz - gradfx1'*diff + (0.5/gam2)*(diff'*diff);
gradFBEx1_gam2 = diff/gam2 - f.Q*diff;

gradfx2 = f.Q*x2 + f.q;
fx2 = 0.5*(x2'*(gradfx2 + f.q));
y = x2 - gam1*gradfx2;
z = proxg(y, gam1);
diff = x2-z;
FBEx2_gam1 = fx2 + gz - gradfx2'*diff + (0.5/gam1)*(diff'*diff);
gradFBEx2_gam1 = diff/gam1 - f.Q*diff;

y = x2 - gam2*gradfx2;
z = proxg(y, gam2);
diff = x2-z;
FBEx2_gam2 = fx2 + gz - gradfx2'*diff + (0.5/gam2)*(diff'*diff);
gradFBEx2_gam2 = diff/gam2 - f.Q*diff;

%% lasso

A = [ 4,    -5,     0,    -3,     1;
		-4,     2,     3,     8,    -1;
	   -11,    -5,     6,    -6,     4;
		 0,     7,   -10,    -1,    -7;
		14,     4,     6,    -6,    -3;
		-2,     5,    -2,     3,   -11;
		-2,    -5,    -8,     2,     1;
		 0,    -7,     5,     1,    -2;
		 0,    -2,    -9,    -2,    -5;
		-5,    -6,    -3,   -11,     4]';
b = [-1, -4, 6, -2, -3]';
gam1 = 0.0017;
gam2 = 0.003;
x1 = [-0.14, -0.24, -0.15, 0.03, 0.03, 0.04, -0.02, 0.01, -0.05, 0.08]';
x2 = [-0.14, 0.04,  -0.09, -0.04,0.05, 0.10, -0.12, 0.12, 0.06, -0.01]';
g = l1Norm(1);
proxg = g.makeprox();

resx1 = A*x1 - b;
gradfx1 = A'*resx1;
fx1 = 0.5*(resx1'*resx1);
y = x1 - gam1*gradfx1;
[z, gz] = proxg(y, gam1);
diff = x1-z;
FBEx1_gam1 = fx1 + gz - gradfx1'*diff + (0.5/gam1)*(diff'*diff);
gradFBEx1_gam1 = diff/gam1 - A'*(A*diff);

y = x1 - gam2*gradfx1;
[z, gz] = proxg(y, gam2);
diff = x1-z;
FBEx1_gam2 = fx1 + gz - gradfx1'*diff + (0.5/gam2)*(diff'*diff);
gradFBEx1_gam2 = diff/gam2 - A'*(A*diff);

% resx2 = A*x2 - b;
% gradfx2 = A'*resx2;
% fx2 = 0.5*(resx2'*resx2);
% y = x2 - gam1*gradfx2;
% [z, gz] = proxg(y, gam1);
% diff = x2-z;
% FBEx2_gam1 = fx2 + gz - gradfx2'*diff + (0.5/gam1)*(diff'*diff);
% gradFBEx2_gam1 = diff/gam1 - A'*(A*diff);
% 
% y = x2 - gam2*gradfx2;
% [z, gz] = proxg(y, gam2);
% diff = x2-z;
% FBEx2_gam2 = fx2 + gz - gradfx2'*diff + (0.5/gam2)*(diff'*diff);
% gradFBEx2_gam2 = diff/gam2 - A'*(A*diff);

%% sparselogreg

% A = [1, -2, 3, -4, 5; 2, -1, 0, -1, 3; -1, 0, 4, -3, 2; -1, -1, -1, 1, 3];
% b = [1, -1, 1, -1]';
% gam1 = 0.1;
% gam2 = 0.05;
% x1 = zeros(5,1);
% f = logLoss(1);
% g = l1Norm(1);
% callf = f.makef();
% proxg = g.makeprox();
% 
% resx1 = A*x1-b;
% [fx1, gradfresx1] = callf(resx1);
% gradfx1 = A'*gradfresx1;
% y = x1 - gam1*gradfx1;
% [z, gz] = proxg(y, gam1);
% diff = x1-z;
% FBEx1_gam1 = fx1 + gz - gradfx1'*diff + (0.5/gam1)*(diff'*diff);
% 
% resx1 = A*x1-b;
% [fx1, gradfresx1] = callf(resx1);
% gradfx1 = A'*gradfresx1;
% y = x1 - gam2*gradfx1;
% [z, gz] = proxg(y, gam2);
% diff = x1-z;
% FBEx1_gam2 = fx1 + gz - gradfx1'*diff + (0.5/gam2)*(diff'*diff);
% 
