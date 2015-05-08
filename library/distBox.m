function obj = distBox(lb,ub,weights)
% Proximal mapping for (weighted) distance from a box [lb,ub]
if nargin<3 || isempty(weights)
    weights = 1;
end
if nargin<2 || isempty(ub)
    ub = +inf;
end
if nargin<1 || isempty(lb)
    lb = -inf;
end

obj.makeprox = @() @(x, gam) call_distBox_prox(x, gam, lb, ub, weights);
end

function [prox, val] = call_distBox_prox(x, gam, lb, ub, weights)
% Proximal mapping of function g(x) = -weights.*min{0,z}
% project on the box
proj = max(min(x,ub),lb);
n = length(x);
if isscalar(weights)
    weights = weights*ones(n,1);
end

if isscalar(lb)
    lb = lb*ones(n,1);
end

if isscalar(ub)
    ub = ub*ones(n,1);
end

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
        prox(iLarge,1) = x(iLarge,1)+gam(iLarge,1).*(diff(iLarge,1)./dist(iLarge,1));
        val = sum(weights(iLarge,1).*abs(max(min(prox(iLarge,1),ub(iLarge,1)),lb(iLarge,1))-prox(iLarge,1)));
    else
        val = 0;
    end
end
end
