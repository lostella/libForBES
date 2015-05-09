function obj = proj_over_affine(p, A, b)
% calculates conjugate function for a quadratic over and affine subspace,
% i.e, solve minimize 0.5*||x-p||^2-y'*x subject to Ax=b
% 
    obj.makefconj = @() make_proj_over_affine_conj(p, A, b);
end


function fc = make_proj_over_affine_conj(p, A, b)
    LD = ldlchol(A,1e-12);
    fc = @(y) eval_proj_over_affine_conj(y, LD, p, A, b);
end

function [val,grad] = eval_proj_over_affine_conj(y, LD, p, A, b)
    yp = y+p;
    grad = yp-A'*ldlsolve(LD,A*yp-b);
    val = y'*grad-0.5*norm(grad-p)^2;
end
