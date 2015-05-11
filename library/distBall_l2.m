function obj = distBall_l2(rho,c,weight)
% Proximal mapping of w*dist(x,C) where C is the l2 ball ||x-c||<=rho
if nargin<3 || isempty(weight)
    weight = 1;
end
if nargin<2 || isempty(c)
    c = 0;
end
if nargin<1 || isempty(rho)
    rho = 1;
end

obj.makeprox = @() @(x, gam) call_distBall_l2_prox(x, gam, rho, c, weight);
end

function [prox, val] = call_distBall_l2_prox(x, gam, rho, c, weight)
gam = weight*gam;
xc = x-c;
nxc = norm(xc);
if nxc <= rho
    prox = x;
    val = 0;
else
    scale = rho/nxc;
    if nxc > gam/(1-scale);
        prox = x-(gam/nxc)*xc;
        xc = prox-c;
        nxc = norm(xc);
        if nxc <= rho
            val = 0;
        else
            val = weight*(1-rho/nxc)*nxc;
        end
    else
        prox = c + scale*xc;
        val = 0;
    end
end
end