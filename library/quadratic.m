function obj = quadratic(Q, q)
% constructs a quadratic 0.5*x'*Qx+q'*x and its conjugate
    obj.makef = @() @(x) call_quadratic(Q, q, x);
    obj.makefconj = @() make_quadratic_conj(Q, q);
end

function [v, g] = call_quadratic(Q, q, x)
    g = Q*x+q;
    v = 0.5*(g+q)'*x;   
end

function fc = make_quadratic_conj(Q, q)
    if issparse(Q)
        [L,flag,p] = chol(Q,'lower','vector');
        if flag~=0
            error('Q is not positive definite')
        end
        fc = @(y) eval_quadratic_sparse_conj(L, p, q, y);
    else
        [L,flag] = chol(Q,'lower');
        if flag~=0
            error('Q is not positive definite')
        end
        fc = @(y) eval_quadratic_dense_conj(L, q, y);
    end
end

function [v, g] = eval_quadratic_dense_conj(L, q, y)
    g = L'\(L\(y-q));
    v = (y-q)'*g;
end

function [v, g] = eval_quadratic_sparse_conj(L, p, q, y)
    rhs = y-q;
    g(p,1) = L'\(L\rhs(p));
    v = (y-q)'*g;
end
