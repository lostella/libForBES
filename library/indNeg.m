function obj = indNeg()
    obj.makeprox = @() @(x, gam) call_indNeg_prox(x);
end

function [prox, val] = call_indNeg_prox(x)
    prox = min(0, x);
    val = 0;
end
