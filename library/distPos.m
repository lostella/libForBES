function obj = distPos(weights,lb)
% Proximal mapping for (weighted) distance from a box [lb,ub]
if nargin<2 || isempty(lb)
    lb = 0;
    if nargin<1 || isempty(weights)
        weights = 1;
    end
end

obj.makeprox = @() @(x, gam) call_distPos_prox(x, gam, lb, weights);
end

function [prox, val] = call_distPos_prox(x, gam, lb, weights)
% Proximal mapping of function g(x) = -weights.*min{0,z}
% project on the box
proj = max(x,lb);
wInf = (weights == inf);
if all(wInf)
    prox = proj;
    val = 0;
else
    diff = proj-x;
    prox = proj;
    gam = gam*weights;
    dist = weights.*abs(diff);
    iLarge = (dist > gam) & ~wInf;
    if any(iLarge)
        prox(iLarge,1) = x(iLarge,1)+gam*(diff(iLarge,1)./dist(iLarge,1));
        if isscalar(lb)
            val = sum(weights(iLarge,1).*abs(max(prox(iLarge,1),lb)-prox(iLarge,1)));
        else            
            val = sum(weights(iLarge,1).*abs(max(prox(iLarge,1),lb(iLarge,1))-prox(iLarge,1)));
        end
    else
        val = 0;
    end
end
end
