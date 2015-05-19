%INDNEG Indicator function of the negative orthant.
%
%   INDNEG(ub) builds the function
%       
%       g(x) = 0    if x <= ub
%            = +inf otherwise
%
%   where ub is either a scalar or a vector of the same size of x. If
%   argument ub is not given, then ub = 0.
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

function obj = indNeg(ub)
    if nargin < 1 || isempty(ub)
        ub = 0;
    end
    obj.makeprox = @() @(x, gam) call_indNeg_prox(x, ub);
end

function [prox, val] = call_indNeg_prox(x, ub)
    prox = min(ub, x);
    val = 0;
end
