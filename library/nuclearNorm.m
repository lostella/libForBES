function obj = nuclearNorm(m, n, lam)
    if nargin < 2
        error('you must provide the number of rows and columns, m and n, as arguments');
    end
    if nargin < 3
        lam = 1;
    end
    obj = @() @(x, gam) call_nuclearNorm_prox(x, gam, m, n, lam);
end

function [y, v] = call_nuclearNorm_prox(x, gam, m, n, lam)
    [U, S, V] = svd(reshape(x, m, n), 'econ');
    diagS1 = max(0, diag(S)-lam*gam);
    S1 = diag(sparse(diagS1));
    y = reshape(U*(S1*V'), m*n, 1);
    if nargout >= 2
        v = lam*sum(diagS1);
    end
end
