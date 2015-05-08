function obj = indBall_l2(rho, c)
% Projection on l2 ball ||x-c||<=rho
if nargin<2 || isempty(c)
    c = 0;
    if nargin<1 || isempty(rho)
        rho = 1;
    end
end
obj.makeprox = @() @(x, gam) call_indBall_l2_prox(x, rho, c);
end

function [prox, val] = call_indBall_l2_prox(x, rho, c)
xc = x - c;
nxc = norm(xc);
if nxc <= rho
    prox = x;
else
    prox = c + (rho/nxc)*xc;
end
val = 0;
end