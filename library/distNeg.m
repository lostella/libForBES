function obj = distNeg(weights,ub)
% Proximal mapping for (weighted) distance from a box [lb,ub]
if nargin<2 || isempty(ub)
    ub = 0;
    if nargin<1 || isempty(weights)
        weights = 1;
    end
end

obj.makeprox = @() @(x, gam) call_distNeg_prox(x, gam, ub, weights);
end

function [prox, val] = call_distNeg_prox(x, gam, ub, weights)
% Proximal mapping of function g(x) = -weights.*min{0,z}
% project on the box
proj = min(x,ub);
wInf = (weights == inf);
if all(wInf)
    prox = proj;
    val = 0;
else
    diff = proj-x;
    prox = proj;
    dist = weights.*abs(diff);
    iLarge = (dist > gam) & ~wInf;
    if any(iLarge)
        prox(iLarge,1) = x(iLarge,1)+gam*(diff(iLarge,1)./dist(iLarge,1));
        if isscalar(ub)
            val = sum(weights(iLarge,1).*abs(min(prox(iLarge,1),ub)-prox(iLarge,1)));
        else            
            val = sum(weights(iLarge,1).*abs(min(prox(iLarge,1),ub(iLarge,1))-prox(iLarge,1)));
        end
    else
        val = 0;
    end
end
end
