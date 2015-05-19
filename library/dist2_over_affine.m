%DIST2_OVER_AFFINE Allocates the squared distance function over an affine subspace.
%
%   DIST2_OVER_AFFINE(p, A, b) returns the function
%       
%       f(x) = 0.5*||x-p||^2 subject to A*x = b
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

function obj = dist2_over_affine(p, A, b)
    obj.makefconj = @() make_dist2_over_affine_conj(p, A, b);
end


function fc = make_dist2_over_affine_conj(p, A, b)
    LD = ldlchol(A,1e-12);
    fc = @(y) call_dist2_over_affine_conj(y, LD, p, A, b);
end

function [val,grad] = call_dist2_over_affine_conj(y, LD, p, A, b)
    yp = y+p;
    grad = yp-A'*ldlsolve(LD,A*yp-b);
    val = y'*grad-0.5*norm(grad-p)^2;
end
