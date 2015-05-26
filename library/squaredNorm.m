%SQUAREDNORM Allocates the squared norm function.
%
%   SQUAREDNORM(lam, p) builds the function
%       
%       f(x) = (lam/2)*||x-p||^2
%   
%   If the arguments are omitted, it is assumed that lam = 1, p = 0.
%
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

function obj = squaredNorm(lam, p)
    if nargin < 1
        lam = 1;
    end
    obj.isQuadratic = 1;
    obj.isConjQuadratic = 1;
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
    val = 0.5*(y'*(grad + p));
end
