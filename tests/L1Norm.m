function [prox, g] = L1Norm(x, gam, lam)
% prox(x) and g(prox(x)) for function g(x) = lam*||x||_1
    uz = max(0, abs(x)-gam*lam);
    if nargout >= 2
        g = lam*sum(uz);
    end
    prox = sign(x).*uz;
end
