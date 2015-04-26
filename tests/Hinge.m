function [prox, g] = Hinge(x, gam, b)
    bx = b.*x; ind = bx < 1;
    prox(ind,1) = b(ind).*min(bx(ind)+gam,1);
    prox(~ind,1) = x(~ind);
    g = sum(max(0,1-b.*prox));
end