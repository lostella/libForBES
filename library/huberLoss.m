function obj = huberLoss(del)
    obj.makef = @() @(x) call_huberLoss_f(x, del);
%     obj.L = 1/del;
end

function [val, grad] = call_huberLoss_f(x, del)
    absx = abs(x);
    small = absx <= del;
    large = ~small;
    sqx = (0.5/del)*(x(small).^2);
    linx = absx(large)-0.5*del;
    val = sum(sqx)+sum(linx);
    if nargout >= 2
        grad = zeros(length(x),1);
        grad(small) = x(small)/del;
        grad(large) = sign(x(large));
    end
end
