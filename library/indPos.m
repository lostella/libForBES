function obj = indPos()
    obj.makeprox = @() @(x, gam) call_indPos_prox(x);
end

function [prox, val] = call_indPos_prox(x)
    prox = max(0, x);
    val = 0;
end
