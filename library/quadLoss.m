%QUADLOSS Allocates the squared norm function.
%
%   QUADLOSS(w, p) builds the function
%       
%       f(x) = (1/2)*||x-p||_Q^2
%   
%   If the arguments are omitted, it is assumed that w = 1, p = 0.
%   If w is a positive scalar then Q = w*Id; if w is a nonnegative vector
%   then Q = diag(w).
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

function obj = quadLoss(w, p)
    if nargin < 1, w = 1; end
    if nargin < 2, p = 0;
    else obj.q = -p; end
    if any(w < 0)
        error('second argument should be a nonnegative');
    end
    if isscalar(w)
        obj.Q = w;
    elseif isvector(w)
        n = length(w);
        obj.Q = spdiags(w,0,n,n);
    end
    obj.isQuadratic = 1;
    obj.isConjQuadratic = 1;
    if all(w > 0)
        obj.makefconj = @() @(x) call_squaredWeightedDistance_conj(x, w, p);
    end
end

function [val, grad] = call_squaredWeightedDistance_conj(y, w, p)
    grad = p + (y./w);
    val = 0.5*(y'*(grad + p));
end
