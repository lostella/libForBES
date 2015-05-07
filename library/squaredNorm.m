function obj = squaredNorm(lam)
    if nargin < 1
        lam = 1;
    end
    obj.makef = @() @(x) call_squaredNorm(x, lam);
    obj.makefconj = @() @(x) call_squaredNorm(x, 1/lam);
end

function [val, grad] = call_squaredNorm(x, lam)
    grad = lam*x;
    val = 0.5*lam*(grad'*x);
end
