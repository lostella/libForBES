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

function obj = distNeg(weights,ub)
    % Proximal mapping for (weighted) distance from a box [lb,ub]
    if nargin<2 || isempty(ub)
        ub = 0;
        if nargin<1 || isempty(weights)
            weights = 1;
        end
    end
    obj.makeprox = @(gam0) @(x, gam) call_distNeg_prox(x, gam, ub, weights);
end

function [prox, val] = call_distNeg_prox(x, gam, ub, weights)
    mu = gam*weights;
    prox = min(max(x-mu,ub),x);
    if nargin>1
        finw = ~isinf(weights);
        val = sum(weights(finw).*max(prox(finw)-ub(finw),0));
    end
end
