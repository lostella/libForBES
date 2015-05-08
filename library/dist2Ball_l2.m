function obj = dist2Ball_l2(rho,c,weight)
% Function value and gradient of (w/2)*dist^2(x,C) where C is the l2 ball ||x-c||<=rho
if nargin<3 || isempty(weight)
    weight = 1;
    if nargin<2 || isempty(c)
        c = 0;
        if nargin<1 || isempty(rho)
            rho = 1;
        end
    end
end
obj.makef = @() @(x) call_dist2Ball_l2_f(x, rho, c, weight);
end

function [val, grad] = call_dist2Ball_l2_f(x, rho, c, weight)
xc = x-c;
nxc = norm(xc);
if nxc <= rho
    proj = x;
else
    proj = c + (rho/nxc)*xc;
end
diff = x-proj;
val = (0.5*w)*(diff'*diff);
grad = weight*diff;
end