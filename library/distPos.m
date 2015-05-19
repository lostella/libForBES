% Copyright (C) 2015, Lorenzo Stella and Panagiotis Patrinos
%
% This file is part of ForBES.
% 
% ForBES is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% ForBES is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU Lesser General Public License for more details.
% 
% You should have received a copy of the GNU Lesser General Public License
% along with ForBES. If not, see <http://www.gnu.org/licenses/>.

function obj = distPos(weights,lb)
% Proximal mapping for (weighted) distance from a box [lb,ub]
if nargin<2 || isempty(lb)
    lb = 0;
end
if nargin<1 || isempty(weights)
    weights = 1;
end

obj.makeprox = @() @(x, gam) call_distPos_prox(x, gam, lb, weights);
end

function [prox, val] = call_distPos_prox(x, gam, lb, weights)
% Proximal mapping of function g(x) = -weights.*min{0,z}
% project on the box
proj = max(x,lb);
n = length(x);
if isscalar(weights)
    weights = weights*ones(n,1);
end

if isscalar(lb)
    lb = lb*ones(n,1);
end


wInf = (weights == inf);
if all(wInf)
    prox = proj;
    val = 0;
else
    diff = proj-x;
    prox = proj;
    gam = gam*weights;
    dist = weights.*abs(diff);
    iLarge = (dist > gam) & ~wInf;
    if any(iLarge)
        prox(iLarge,1) = x(iLarge,1)+gam(iLarge,1).*(diff(iLarge,1)./dist(iLarge,1));
        val = sum(weights(iLarge,1).*abs(max(prox(iLarge,1),lb(iLarge,1))-prox(iLarge,1)));
    else
        val = 0;
    end
end
end
