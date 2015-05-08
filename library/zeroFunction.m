function obj = zeroFunction()
    obj.makeprox = @() @(x, gam) call_zeroFunction_prox(x);
end

function [prox, val] = call_zeroFunction_prox(x)
    prox = x;
    val = 0;
end