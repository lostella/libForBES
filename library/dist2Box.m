function obj = dist2Box(lb, ub, weights)
    % Function value and gradient of (w/2)*dist^2(x,C) where C is the box [lb,ub]
    if nargin<3 || isempty(weights)
        weights = 1;
    end
    if nargin<2 || isempty(ub)
        ub = +inf;
    end
    if nargin<1 || isempty(lb)
        lb = -inf;
    end

    obj.makef = @() @(x) call_dist2Box_f(x, lb ,ub, weights);
    obj.L = max(weights);
end

function [val, grad] = call_dist2Box_f(x,lb ,ub, weights)
    proj = max(min(x,ub),lb);
    diff = x - proj;
    grad = weights.*diff;
    val  = 0.5*sum(weights.*(diff.*diff));
end