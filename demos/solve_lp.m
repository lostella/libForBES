function [x, y, fval, out] = solve_lp(c,A,b)
    [m,n] = size(A);
    K = [sparse(n, n), A', speye(n); A, sparse(m, m+n); c', -b', sparse(1, n)];
    d = [c; b; 0];
    f = quadLossOverAffine(sparse(2*n+m,1),K,d);
    g = indPos([repmat(0, n, 1); repmat(-inf, m, 1); repmat(0, n, 1)]);
    constr = {1, -1, zeros(2*n+m,1)};
    out = forbes(f, g, zeros(2*n+m,1), [], constr);
    x = out.x1(1:n);
    y = out.x1(n+1:n+m);
    fval = c'*x;
end
