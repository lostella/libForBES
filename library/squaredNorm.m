function obj = squaredNorm(lam,p)
% constructs (lam/2)||x-p||^2 and its conjugate
% arguments lam p are optional. if omitted then
% lam = 1 and p=0;
    if nargin < 1
        lam = 1;
    end
    if nargin==2 && ~isempty(p)
        obj.makef = @() @(x) call_squaredDist(x, lam, p);
        obj.makefconj = @() @(x) call_squaredDist_conj(x, 1/lam, p);
    else
        obj.makef = @() @(x) call_squaredNorm(x, lam);
        obj.makefconj = @() @(x) call_squaredNorm(x, 1/lam);
    end
end

function [val, grad] = call_squaredNorm(x, lam)
    grad = lam*x;
    val = 0.5*(grad'*x);
end

function [val, grad] = call_squaredDist(x, lam, p)
    xp = x - p;
    grad = lam*xp;
    val = 0.5*(grad'*xp);
end

function [val, grad] = call_squaredDist_conj(y, lam, p)
    grad = p + lam*y;
    val =0.5*(y'*(grad + p));
end