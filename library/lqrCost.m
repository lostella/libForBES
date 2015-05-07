function obj = lqrCost(x0, A, B, LRs, Ks, Ms, Ls, N)
    obj.makefconj = @() make_lqrCost_fconj(x0, A, B, LRs, Ks, Ms, Ls, N);
end

function op = make_lqrCost_fconj(x0, A, B, LRs, Ks, Ms, Ls, N)
    [n, m] = size(B);
    op = @(w) RiccatiSolve(w, x0, A, B, LRs, Ks, Ms, Ls, int32(n), int32(m), int32(N));
end
