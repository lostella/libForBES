%INDFREEPOS Partial indicator function of the positive orthant.
%
%   INDFREEPOS(iPos) builds the function
%       
%       g(x) = 0    if x_i >= 0 for all i in iPos
%            = +inf otherwise
%
%   where iPos is a vector containing integers.
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

function obj = indFreePos(iPos)
    obj.makeprox = @() @(x, gam) call_indFreePos_prox(x, iPos);
end

function [prox, val] = call_indFreePos_prox(x, iPos)
    prox = x;
    prox(iPos,1)=max(0, x(iPos,1));
    val = 0;
end
