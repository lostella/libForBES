%LOGLOGISTIC Allocates the log-logistic function.
%
%   LOGLOGISTIC(mu) builds the log-logistic function
%       
%       f(x) = mu*(sum_i log(1+exp(-x_i)))
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

function obj = logLogistic(mu)
    obj.makef = @() @(x) call_logLogistic_f(x, mu);
    obj.L = mu; % Lipschitz constant of the gradient of f
end

function [val, grad] = call_logLogistic_f(x, mu)
% value and gradient of f(x) = mu*sum(log(1+exp(-x)))
    px = 1./(1+exp(-x));
    val = -sum(log(px))*mu;
    if nargout >= 2
        grad = (px-1)*mu;
    end
end
