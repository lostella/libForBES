%HINGELOSS Allocates the hinge loss function.
%
%   HINGELOSS(mu, b) builds the function
%       
%       g(x) = mu*sum(max(0, 1-b.*x))
%

function obj = hingeLoss(mu, b)
    %
    % Only the proximal mapping is available for this function
    %
    obj.makeprox = @() @(x, gam) call_hingeLoss_prox(x, gam, mu, b);
end

function [prox, g] = call_hingeLoss_prox(x, gam, mu, b)
    bx = b.*x; ind = bx < 1;
    prox(ind,1) = b(ind).*min(bx(ind)+gam*mu,1);
    prox(~ind,1) = x(~ind);
    g = sum(max(0,1-b.*prox));
end

