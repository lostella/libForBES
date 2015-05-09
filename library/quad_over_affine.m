function obj = quad_over_affine(Q, q, A, b)
% calculates conjugate function for a quadratic over and affine subspace,
% i.e, solve minimize 0.5*x'*Q*x+q'*x-y'*x subject to Ax=b
% 
    obj.makefconj = @() make_quad_over_affine_conj(Q, q, A, b);
end


function fc = make_quad_over_affine_conj(Q, q, A, b)
    m = size(A,1);
    [L,D,P] = ldl([Q A';A sparse(m,m)]);
    fc = @(y) eval_quad_over_affine_conj(y, L, D, P, Q, q, b);
end

function [val, grad] = eval_quad_over_affine_conj(y, L, D, P, Q, q, b)
    grad = P*(L'\(D\(L\(P'*[y-q; b]))));
    grad = grad(1:length(q),1);
    val = -(0.5*grad'*(Q*grad)+(q-y)'*grad);
end
