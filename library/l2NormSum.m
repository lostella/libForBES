function obj = l2NormSum(m)
    if nargin < 1
        m = 1;
    end
    obj = @() @(x, gam) call_l2NormSum_prox(x, gam, m);
end

function [z, v] = call_l2NormSum_prox(x, gam, m)
% prox(x) and g(prox(x)) for function g(x) being the sum of the Euclidean
% norm in R^m of the blocks of size m of x
    n = length(x);
    nb = n/m;
    x_resh = reshape(x, m, nb);
    modx = sqrt(sum(x_resh.*x_resh, 1));
    newmodx = max(0, modx-gam);
    scal_block = newmodx./max(gam, modx);
    scal_resh = repmat(scal_block, m, 1);
    scal = reshape(scal_resh, n, 1);
    z = scal.*x;
    if nargout >= 2
        v = sum(newmodx);
    end
end
