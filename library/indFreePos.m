function obj = indFreePos(n,m)
% Proximal mapping of indicator of 
% C={x_i, 1 = 1,...,n,n+1,...n+m | x_i, i = 1,...,n free, x_i >= 0, i=n+1,...,n+m}
    obj.makeprox = @() @(x, gam) call_indPos_prox(x,n,m);
end

function [prox, val] = call_indPos_prox(x,n,m)
    prox = [x(1:n,1);max(0, x(n+1:n+m,1))];
    val = 0;
end
