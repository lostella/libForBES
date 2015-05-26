%FORBES Solver for nonsmooth convex optimization problems.
% 
%   We assume the problem at hand has the form
%
%           minimize    f1(C1 x1 - d1) + f2(C2 x2 - d2) + g(z),
%           subject to  A1 x1 + A2 x2 + Bz = b
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
%   If the constraint is not specified, it is assumed to be x1 = x2 = z, in
%   which case the problem takes the unconstrained form
%
%       (1) minimize    f1(C1 x - d1) + f2(C2 x - d2) + g(x)
%
%   and we assume that both f1 and f2 are smooth, and f1 is quadratic.
%   If A, B and b are specified then C = Id and d = 0, and the problem has
%   the form
%
%       (2) minimize    f1(x1) + f2(x2) + g(z)
%           subject to  A1 x1 + A2 x2 + Bz = b
%
%   And we assume that f1 and f2 are strongly convex, f1 being quadratic
%   plus (at most) the indicator of an affine subspace (so that its
%   conjugate is quadratic).
%
%   The problem is defined through the following attributes of the
%   structure prob provided as first argument to FORBES.
%
%       prob.x0: The starting point for the algorithm.
% 
%       prob.f1, prob.f2, prob.g: functions f and g in the cost. See section
%           "Functions" below for details on how to select this.
%           (One of f1 and f2 is required, g is required).
%
%       prob.C1, prob.d1: Matrix (or function handle) and vector with which
%           f1 is composed. (Both are optional; default: C1 = Id, d1 = 0).
%
%       prob.C2, prob.d2: Matrix (or function handle) and vector with which
%           f2 is composed. (Both are optional; default: C2 = Id, d2 = 0).
% 
%       prob.C1t, prob.C2t: Function computing the adjoint of C1 and C2.
%           (Required if C1 and C2 are specified as function handles).
% 
%       prob.A1, prob.A2: matrix (or procedures computing matvecs) defining
%           the constraint. (Optional).
% 
%       prob.A1t, prob.A2t: procedures computing the adjoint of A1 and A2.
%           (Required only if A1 or A2 are specified as function handles).
%
%       prob.B: matrix B in the equality constraint.
%           (Optional).
% 
%       prob.b: right hand side vector in the constraint.
%           (Optional).
%
%   Functions
%   ---------
%
%   Functions f and g in the cost can be selected in a library of functions
%   available in the "library" directory inside of FORBES directory. All
%   these functions return a structure containing all that is needed by the
%   solver, and may accept parameters defining the function. For example
%
%       prob.f = logLogistic(mu)
%       prob.g = l1Norm(mu)
%
%   puts in prob.f the log-logistic loss function
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
%       opt.tol: Tolerance on the optimality condition. (Default: 1e-8).
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
% Copyright (C) 2015, Lorenzo Stella and Panagiotis Patrinos
%
% This file is part of ForBES.
% 
% ForBES is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% ForBES is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU Lesser General Public License for more details.
% 
% You should have received a copy of the GNU Lesser General Public License
% along with ForBES. If not, see <http://www.gnu.org/licenses/>.

function out = forbes(prob, opt)
    prob.id = IdentifyProblem(prob);
    switch prob.id
        case 1
            if nargin > 1
                if isfield(opt, 'method') && strcmp(opt.method, 'fbs'), out = fbs(prob, opt);
                else out = minfbe(prob, opt); end
            else
                out = minfbe(prob);
            end
        case 2
            if nargin > 1
                if isfield(opt, 'method') && strcmp(opt.method, 'fbs'), out = amm(prob, opt);
                else out = miname(prob, opt); end
            else
                out = miname(prob);
            end
    end
end

function id = IdentifyProblem(prob)
    % simply check for the presence of linear equality constraints.
    if any(isfield(prob, {'A', 'At', 'B', 'b'}))
        flagEq = true;
    else
        flagEq = false;
    end
    if flagEq
        id = 2;
    else
        id = 1;
    end
end
