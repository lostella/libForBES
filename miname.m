%MINAME Solver for equality constrained convex problems
% 
%   We assume the problem at hand has the form
% 
%       (1) minimize f1(x1) + f2(x2) + g(z)
%           subject to A1x1 + A2x2 + Bz = c.
%
%   Moreover, f1 and f2 are strongly convex and twice continuously
%   differentiable, with f1 being quadratic the indicator of an affine
%   subspace S (which may be the whole R^n), that is 
%
%       f1(x1) = 1/2 x1'Qx1 + q'x1, if x1 in S
%              = +infinity, otherwise
%
%   OUT = MINAME(PROB) solves problem (1) specified by the structure PROB,
%   using the default options, and putting the results of the optimization
%   process in OUT.
%
%   OUT = MINAME(PROB, OPT) like the previous, but uses the options specified
%   in the OPT structure instead of the defaults.
%
%   The problem structure PROB
%   --------------------------
%
%   The algorithm requires in the PROB argument, defining the problem, the
%   following fields:
%
%       PROB.x1step: a procedure solving, given w, the problem
% 
%               x1 = argmin{f1(x) - w'x}.
% 
%           This corresponds to computing the gradient of the conjugate of
%           f1 (Optional).
% 
%       PROB.x1step: a procedure solving, given w, the problem
% 
%               x2 = argmin{f2(x) - w'x}.
% 
%           This corresponds to computing the gradient of the conjugate of
%           f2 (Optional).
% 
%       PROB.zstep: procedure minimizing, given y and gamma, the augmented
%       Lagrangian, with parameter gamma, with respect to z. That is
% 
%              z = argmin{g(z) + y'Bz + gamma/2 ||Bz||^2},
% 
%           and g(z) is returned as second output. (Required).
% 
%       PROB.A1, PROB.A2, PROB.B: matrices (or procedures computing matvecs)
%           defining the constraint. (Required if the corresponding term in
%           the cost is specified).
% 
%       PROB.A1T, PROB.A2T, PROB.BT: procedures computing the adjoint of
%           A1, A2, B. (Required only if A1 or A2 or B are specified as
%           function handles).
% 
%       PROB.c: right hand side vector in the constraint. (Required).
%
%   The options structure OPT
%   -------------------------
%
%   In OPT the user can specify the behaviour of the algorithm to be used.
%   The following options can be set:
%
%       OPT.tolOpt: Tolerance on the optimality condition. (Default: 1e-5).
% 
%       OPT.maxit: Maximum number of iterations. (Default: 10*n, with n
%           being the number of variables).
% 
%       OPT.method: Algorithm to use. Can select between:
%           * 'sd' (steepest descent),
%           * 'lbfgs' (limited memory BFGS, default),
%           * 'cg-desc', 'cg-prp', 'cg-dyhs' (various CG algorithms),
%           * 'bb' (Barzilai-Borwein).
% 
%       OPT.variant: Variant of the method to use. Select between:
%           * 'basic', the basic algorithm,
%           * 'global', the global variant (default),
%           * 'fast', the fast variant.
% 
%       OPT.linesearch: Line search strategy to use. Can select between:
%           * 'armijo' (default for 'sd'),
%           * 'nonmonotone-armijo' (default for 'bb'),
%           * 'hager-zhang' (default for the rest),
%           * 'lemarechal',
%           * 'fletcher'.
%
%   The output
%   ----------
%   
%   The OUT structure contains, besides the solution, a series of details
%   concerning the optimization process. The main attributes of OUT are the
%   following:
% 
%       OUT.name: string summarizing the method used.
% 
%       OUT.message: output message.
% 
%       OUT.flag: output flag. Can be:
%           * 0, if the solver converged up to the prescribed tolerance,
%           * 1, if the maximum number of iterations was exceeded,
%           * 2, if the line search failed.
% 
%       OUT.y: computed dual solution found.
% 
%       OUT.x1, OUT.x2, OUT.z: computed primal solution.
% 
%       OUT.iterations: number of iterations taken.
% 
%       OUT.operations: summary of the number of operations performed.
% 
%       OUT.residual: evolution of the residual during the iterates.
% 
%       OUT.ts: timestamp for each of the iterations.
% 
%   See also MINFBE
%
% Authors: Lorenzo Stella (lorenzo.stella -at- imtlucca.it)
%          Panagiotis Patrinos (panagiotis.patrinos -at- imtlucca.it)

function out = miname(prob, opt)
    t0 = tic();
    dualprob = ProcessDualProblem(prob);
    preprocess = toc(t0);
    if nargin > 1
        dualout = minfbe(dualprob, opt);
    else
        dualout = minfbe(dualprob);
    end
    out = GetPrimalOutput(prob, dualout);
    out.preprocess = preprocess + dualout.preprocess;
    if nargin > 1
        out.opt = opt;
    end
    out.dual = dualout;
end
