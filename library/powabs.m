function obj = powabs(p,c)
% Proximal mapping of eta*|x|^p
% see Combettes, Dung, Vu, Dualization of Signal Recovery Problems
% Example 2.15
if nargin<2 || isempty(c)
    c = 1;
    if nargin<1 || isempty(p)
        p = 1;
    end
end
pows = [1;4/3;3/2;2;3;4];
if all((p == pows) == 0)
    error('Prox is not computable')
end
obj.makeprox = @() @(x, gam) call_powabs_prox(x, gam, p, c);
end

function [prox, val] = call_powabs_prox(x, gam, p, c)
gam = c*gam;
switch p
    case 1
        prox = sign(x)*max(abs(x)-gam,0);
    case 4/3
        rho = sqrt(x^2+256*gam^3/729);
        prox = x + (4*gam)/(3*2^(1/3))*(abs(rho-x)^(1/3)-abs(rho+x)^(1/3));
    case 3/2
        prox = x + 9*gam^2*sign(x)*(1-sqrt(1+16*abs(x)/(9*gam^2)))/8;
    case  2
        prox = x/(1+2*gam);
    case  3
        prox = sign(x)*(sqrt(1+12*gam*abs(x))-1)/(6*gam);
    case  4
        rho = sqrt(x^2+1/(27*gam));
        prox = abs((rho + x)/(8*gam))^(1/3)-abs((rho - x)/(8*gam))^(1/3);
end
val = c*abs(prox)^p;
end