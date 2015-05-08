function obj = dist2Neg(weights,ub)
% Function value and gradient of (w/2)*dist^2(x,C) where C is the box [-infty,ub]
if nargin<2 || isempty(ub)
    ub = 0;
    if nargin<1 || isempty(weights)
        weights = 1;
    end
end
obj.makef = @() @(x) call_dist2Neg_f(x, ub, weights);
end

function [val, grad] = call_dist2Neg_f(x, ub, weights)
proj = min(x,ub);
diff = x - proj;
grad = weights.*diff;
val  = 0.5*sum(weights.*(diff.*diff));
end