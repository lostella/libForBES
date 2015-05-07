function obj = softPositiveOrthant(weights)
    obj.makeprox = @() @(x, gam) call_softPositiveOrthant_prox(x, gam, weights);
end

function [prox, val] = call_softPositiveOrthant_prox(x, gam, weights)
% Proximal mapping of function g(x) = -weights.*min{0,z}
    prox = x;
    neg = prox < 0;
    infw = weights == +inf;
    if isscalar(weights)
        if infw
            prox(neg) = 0;
            val = 0;
        else
            prox(neg) = min(0, x(neg)+gam*weights);
            val = -sum(weights*prox(neg));
        end
    else
        tozero = (neg & infw);
        toshrink = (neg & (~infw));
        prox(tozero) = 0;
        prox(toshrink) = min(0, x(toshrink)+gam*weights(toshrink));
        val = -sum(weights(toshrink).*prox(toshrink));
    end
end
