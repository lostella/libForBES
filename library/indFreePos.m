function obj = indFreePos(iPos)
% Proximal mapping of indicator of 
% C={x_i, 1 = 1,...,n,n+1,...n+m | x_i >= 0, i\in iPos}
    obj.makeprox = @() @(x, gam) call_indPos_prox(x,iPos);
end

function [prox, val] = call_indPos_prox(x,iPos)
    prox = x;
    prox(iPos,1)=max(0, x(iPos,1));
    val = 0;
end
