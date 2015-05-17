%ZEROFUNCTION Allocates the function identically equal to zero.
%
%   ZEROFUNCTION(n) builds the function
%       
%       f(x) = 0
%
%   where x is a vector of length n.
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

function obj = zeroFunction()
    obj.makeprox = @() @(x, gam) call_zeroFunction_prox(x);
end

function [prox, val] = call_zeroFunction_prox(x)
    prox = x;
    val = 0;
end