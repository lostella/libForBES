function obj = hingeLoss(mu, b)
    obj.makeprox = @() @(x, gam) call_hingeLoss_prox(x, gam, mu, b);
end

function [prox, g] = call_hingeLoss_prox(x, gam, mu, b)
    bx = b.*x; ind = bx < 1;
    prox(ind,1) = b(ind).*min(bx(ind)+gam*mu,1);
    prox(~ind,1) = x(~ind);
    g = sum(max(0,1-b.*prox));
end
