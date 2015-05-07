function obj = indBox(lower, upper)
    obj.makeprox = @() @(x, gam) call_indBox_prox(x, lower, upper);
end

function [prox, val] = call_indBox_prox(x, lower, upper)
    prox = min(upper, max(lower, x));
    val = 0;
end
