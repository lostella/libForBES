function obj = distNeg_new(weights,ub)
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
mu = gam*weights;
prox = min(max(x-mu,ub),x);
if nargin>1
    finw = ~isinf(weights);
    val = sum(weights(finw).*max(prox(finw)-ub(finw),0));
end
end
