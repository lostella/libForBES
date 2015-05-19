%INDPOS Indicator function of the positive orthant.
%
%   INDPOS(lb) builds the function
%       
%       g(x) = 0    if x >= lb
%            = +inf otherwise
%
%   where lb is either a scalar or a vector of the same size of x. If
%   argument lb is not given, then lb = 0.
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

function obj = indPos(lb)
    if nargin < 1 || isempty(lb)
        lb = 0;
    end
    obj.makeprox = @() @(x, gam) call_indPos_prox(x, lb);
end

function [prox, val] = call_indPos_prox(x, lb)
    prox = max(lb, x);
    val = 0;
end
