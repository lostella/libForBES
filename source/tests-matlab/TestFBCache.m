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
z = proxg(y, gam1);
diff = x1-z;
gradFBEx1_gam1 = diff/gam1 - f.Q*diff;

y = x1 - gam2*gradfx1;
z = proxg(y, gam2);
diff = x1-z;
gradFBEx1_gam2 = diff/gam2 - f.Q*diff;

%% sparselogreg

A = [1, -2, 3, -4, 5; 2, -1, 0, -1, 3; -1, 0, 4, -3, 2; -1, -1, -1, 1, 3];
b = [1, -1, 1, -1]';
gam1 = 0.1;
gam2 = 0.05;
x1 = zeros(5,1);
f = logLoss(1);
g = l1Norm(1);
callf = f.makef();
proxg = g.makeprox();

resx1 = A*x1-b;
[fx1, gradfresx1] = callf(resx1);
gradfx1 = A'*gradfresx1;
y = x1 - gam1*gradfx1;
[z, gz] = proxg(y, gam1);
diff = x1-z;
FBEx1_gam1 = fx1 + gz - gradfx1'*diff + (0.5/gam1)*(diff'*diff);

resx1 = A*x1-b;
[fx1, gradfresx1] = callf(resx1);
gradfx1 = A'*gradfresx1;
y = x1 - gam2*gradfx1;
[z, gz] = proxg(y, gam2);
diff = x1-z;
FBEx1_gam2 = fx1 + gz - gradfx1'*diff + (0.5/gam2)*(diff'*diff);

