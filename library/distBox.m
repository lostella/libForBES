function obj = distBox(lb,ub,weights)
% Proximal mapping for (weighted) distance from a box [lb,ub]

obj.makeprox = @() @(x, gam) call_distBox_prox(x, gam, lb, ub, weights);
end

function [prox, val] = call_distBox_prox(x, gam, lb, ub, weights)
mu = gam*weights;
prox = max(x - ub - mu,0) - max(lb - x - mu,0) + min(max(x,lb),ub);
if nargout>1
    finw = ~isinf(weights);    
    val = sum(weights(finw).*abs(prox(finw)-min(max(prox(finw),lb(finw)),ub(finw))));
end
end
