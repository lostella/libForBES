function obj = dist2Pos(weights,lb)
% Function value and gradient of (w/2)*dist^2(x,C) where C is the box [lb,+infty]
if nargin<2 || isempty(lb)
    lb = 0;
    if nargin<1 || isempty(weights)
        weights = 1;
    end
end
obj.makef = @() @(x) call_dist2Pos_f(x, lb, weights);
end

function [val, grad] = call_dist2Pos_f(x, lb, weights)
proj = max(x,lb);
diff = x - proj;
grad = weights.*diff;
val  = 0.5*sum(weights.*(diff.*diff));
end