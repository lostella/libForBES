function [fx, gradfx] = Huber(x, del)
    absx = abs(x);
    small = absx <= del;
    large = ~small;
    sqx = (0.5/del)*(x(small).^2);
    linx = absx(large)-0.5*del;
    fx = sum(sqx)+sum(linx);
    if nargout >= 2
        gradfx = zeros(length(x),1);
        gradfx(small) = x(small)/del;
        gradfx(large) = sign(x(large));
    end
end
