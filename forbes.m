%FORBES Solver for nonsmooth convex optimization problems.
% 
%   We assume the problem at hand has the form
%
%           minimize    f1(C1*x1-d1) + f2(C2*x2-d2) + g(z),
%           subject to  A1*x1 + A2*x2 + B*z = b
%
%   out = FORBES(prob) solves problem (1) specified by the structure prob,
%   using the default options, and store the results in out.
%
%   out = FORBES(prob, opt) like the previous, but uses the options specified
%   in the opt structure instead of the defaults.
%
%   Defining the problem
%   --------------------
%
%   If A1, A2, B and b are specified then C1 = C2 = Id and d1 = d2 = 0,
%   and the problem has the form
%
%       (1) minimize    f1(x1) + f2(x2) + g(z)
%           subject to  A1*x1 + A2*x2 + B*z = b
%
%   In this case we assume that f1 is strongly convex and quadratic plus
%   (at most) the indicator of an affine subspace, while f2 is strongly
%   convex in the interior of its domain. Function g is any closed, proper,
%   convex function.
%
%   In case (1) the problem is defined through the following attributes of
%   structure prob.
%
%       prob.f1, prob.f2, prob.g: functions f1, f2, g in the cost, if present.
%           See section "Functions" below for details on how to select these
%           functions.
% 
%       prob.A1, prob.A2: matrices (or procedures computing matvecs)
%           defining the constraint.
% 
%       prob.A1t, prob.A2t: procedures computing the adjoint of
%           A1, A2, B. (Required only if A1 or A2 or B are specified as
%           function handles).
%
%       prob.B: matrix B in the equality constraint.
% 
%       prob.b: right hand side vector in the constraint.
%
%   If the constraint is not specified, it is assumed to be x1 = x2 = z,
%   in which case the problem takes the unconstrained form
%
%       (2) minimize    f1(C1*x-d1) + f2(C2*x-d2) + g(x)
%
%   and we assume that f1 is convex quadratic, f2 is convex, smooth (has
%   Lipschitz continuous gradient) and twice continuously differentiable,
%   while g is any closed, proper, convex function.
%
%   In case (2) the structure prob describing the problem should contain the
%   following attributes.
%
%       prob.x0: The starting point for the algorithm.
% 
%       prob.Q, prob.q: Hessian and linear parts of the quadratic term f1.
%           Q may be a function handle instead of a matrix. (Required only
%           if the quadratic term is present; default: Q = 0, q = 0).
% 
%       prob.C1, prob.d1: Matrix (or function handle) and vector with which
%           f1 is composed. (Both are optional; default: C1 = Id, d1 = 0).
% 
%       prob.C1t: Function computing the adjoint of C1. (Required if C1 is
%           specified as a function handle).
% 
%       prob.f2: function f2 in the cost, if present. See section "Functions"
%           below for details on how to select this.
% 
%       prob.C2, prob.d2: Matrix (or function handle) and vector with which
%           f2 is composed. (Both are optional; default: C2 = Id, d2 = 0).
% 
%       prob.C2t: Function computing the adjoint of C2. (Required if C2 is
%           specified as a function handle).
% 
%       prob.normC2: The 2-norm of matrix C2. (Optional)
% 
%       prob.g: function g in the cost. See section "Functions" below for
%           details on how to select this. (Required).
%
%   Functions
%   ---------
%
%   Functions f1, f2 and g in the cost can be selected in a library of
%   functions available in the "library" directory inside of FORBES
%   directory. All these functions return a structure containing all that
%   is needed by the solver, and may accept parameters defining the
%   function. For example
%
%       prob.f2 = logLogistic(mu)
%       prob.g = l1Norm(mu)
%
%   puts in prob.f2 the log-logistic loss function
%
%       f(x) = mu*(sum_i log(1+exp(-x_i)))
%
%   while prob.g will contain the l1 regularization term mu*||x||_1.
%   Consider looking into the "library" directory for specific information
%   any of the functions.
%
%   The options structure opt
%   -------------------------
%
%   In opt the user can specify the behaviour of the algorithm to be used.
%   The following options can be set:
%
%       opt.tolOpt: Tolerance on the optimality condition. (Default: 1e-5).
% 
%       opt.maxit: Maximum number of iterations. (Default: 10*n, with n
%           being the number of variables).
% 
%       opt.method: Algorithm to use. Can select between:
%           * 'sd' (steepest descent),
%           * 'lbfgs' (limited memory BFGS, default),
%           * 'cg-desc', 'cg-prp', 'cg-dyhs' (various CG algorithms),
%           * 'bb' (Barzilai-Borwein).
% 
%       opt.variant: Variant of the method to use. Select between:
%           * 'basic', the basic algorithm,
%           * 'global', the global variant (default),
%           * 'fast', the fast variant.
% 
%       opt.linesearch: Line search strategy to use. Can select between:
%           * 'armijo' (default for 'sd'),
%           * 'nonmonotone-armijo' (default for 'bb'),
%           * 'hager-zhang' (default for the rest),
%           * 'lemarechal',
%           * 'fletcher'.
%
% Authors: Lorenzo Stella (lorenzo.stella -at- imtlucca.it)
%          Panagiotis Patrinos (panagiotis.patrinos -at- imtlucca.it)
%

function out = forbes(prob, opt)
    prob = IdentifyProblem(prob);
    switch prob.identified
        case 1
            if nargin > 1
                out = minfbe(prob, opt);
            else
                out = minfbe(prob);
            end
        case 2
            if nargin > 1
                out = miname(prob, opt);
            else
                out = miname(prob);
            end
    end
end

function prob = IdentifyProblem(prob)
    % Check whether we are given equality constraints or not
    if any(isfield(prob, {'A1', 'A1t', 'A2', 'A2t', 'B', 'b'}))
        flagEQ = true;
    else
        flagEQ = false;
    end
    % Check whether we are given affine mappings composed with f1 and f2
    if any(isfield(prob, {'C1', 'C1t', 'd1', 'C2', 'C2t', 'd2'}))
        flagAFF = true;
    else
        flagAFF = false;
    end
    % Check for uncertain situations
    if (flagEQ && flagAFF)
        error('you cannot provide both equality constraints and affine mappings composed with f1, f2');
    end
    if flagEQ
        % If equality constraints are provided, solve the dual
        prob.identified = 2;
    else
        % Otherwise solve the primal
        prob.identified = 1;
    end
end
