%LOGLOGISTIC Allocates the log-logistic function.
%
%   LOGLOGISTIC(mu) builds the log-logistic function
%       
%       f(x) = mu*(sum_i log(1+exp(-x_i)))
%

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
