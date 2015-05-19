%QUAD_OVER_AFFINE Allocates a quadratic function over an affine subspace.
%
%   QUADRATIC(Q, q, A, b) builds the function
%       
%       f(x) = 0.5*(x'*Q*x)+q'*x subject to A*x = b
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

function obj = quad_over_affine(Q, q, A, b)
    obj.makefconj = @() make_quad_over_affine_conj(Q, q, A, b);
end


function fc = make_quad_over_affine_conj(Q, q, A, b)
    m = size(A,1);
    [L,D,P] = ldl([Q A';A sparse(m,m)]);
    fc = @(y) eval_quad_over_affine_conj(y, L, D, P, Q, q, b);
end

function [val, grad] = eval_quad_over_affine_conj(y, L, D, P, Q, q, b)
    grad = P*(L'\(D\(L\(P'*[y-q; b]))));
    grad = grad(1:length(q),1);
    val = -(0.5*grad'*(Q*grad)+(q-y)'*grad);
end
