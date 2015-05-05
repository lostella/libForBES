function [prox, gprox] = SoftNonnegativeIndicator(x, gam, weights)
% Proximal mapping of function g(x) = -weights.*min{0,z}
    prox = x;
    neg = prox < 0;
    infw = weights == +inf;
    if isscalar(weights)
        if infw
            prox(neg) = 0;
            gprox = 0;
        else
            prox(neg) = min(0, x(neg)+gam*weights);
            gprox = -sum(weights*prox(neg));
        end
    else
        tozero = (neg & infw);
        toshrink = (neg & (~infw));
        prox(tozero) = 0;
        prox(toshrink) = min(0, x(toshrink)+gam*weights(toshrink));
        gprox = -sum(weights(toshrink).*prox(toshrink));
    end
end

