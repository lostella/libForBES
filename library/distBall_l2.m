function obj = distBall_l2(rho,c,weight)
% Proximal mapping of w*dist(x,C) where C is the l2 ball ||x-c||<=rho
if nargin<3 || isempty(weight)
    weight = 1;
    if nargin<2 || isempty(c)
        c = 0;
        if nargin<1 || isempty(rho)
            rho = 1;
        end
    end
end
obj.makeprox = @() @(x, gam) call_distBall_l2_prox(x, gam, rho, c, weight);
end

function [prox, val] = call_distBall_l2_prox(x, gam, rho, c, weight)
gam = weight*gam;
xc = x-c;
nxc = norm(xc);
if nxc <= rho
    proj = x;
else
    proj = c + (rho/nxc)*xc;
end
diff = proj - x;
dist = norm(diff);
prox = proj;
val = 0;
if dist > gam
    prox = x+gam*(diff/dist);
    xc = prox-c;
    nxc = norm(xc);
    if nxc <= rho
        val = 0;
    else
        val = weight*(1-rho/nxc)*norm(xc);
    end    
end
end