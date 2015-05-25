%FORBES Solver for nonsmooth convex optimization problems.
% 
%   We assume the problem at hand has the form
%
%           minimize    f(Cx-d) + g(z),
%           subject to  Ax + Bz = b
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
%   If the constraint is not specified, it is assumed to be x = z, in which
%   case the problem takes the unconstrained form
%
%       (1) minimize    f(Cx-d) + g(x)
%
%   If A, B and b are specified then C = Id and d = 0, and the problem has
%   the form
%
%       (2) minimize    f(x) + g(z)
%           subject to  Ax + Bz = b
%   
%   The problem is defined through the following attributes of the
%   structure prob provided as first argument to FORBES.
%
%       prob.x0: The starting point for the algorithm.
% 
%       prob.f, prob.g: functions f and g in the cost. See section "Functions"
%           below for details on how to select this. (Required).
%
%       prob.C, prob.d: Matrix (or function handle) and vector with which
%           f is composed. (Both are optional; default: C = Id, d = 0).
% 
%       prob.Ct: Function computing the adjoint of C.
%           (Required if C is specified as a function handle).
% 
%       prob.A: matrix (or procedures computing matvecs) defining the
%           constraint. (Optional; default: A = Id).
% 
%       prob.A1t: procedures computing the adjoint of A.
%           (Required only if A is specified as a function handle).
%
%       prob.B: matrix B in the equality constraint.
%           (Optional; default: B = -Id).
% 
%       prob.b: right hand side vector in the constraint.
%           (Optional; default: b = 0).
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
    t0 = tic();
    prob.id = IdentifyProblem(prob);
    prob = ProcessProblem(prob);
    preprocess = toc(t0);
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
    out = ProcessOutput(out);
    out.preprocess = preprocess + out.preprocess;
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

function obj = ProcessFunction(obj)
    % fill in all fields describing the function with default values
    % (in case they are missing).
    if ~isfield(obj, 'isQuadratic'), obj.isQuadratic = 0; end
    if ~isfield(obj, 'isConjQuadratic'), obj.isConjQuadratic = 0; end
    if ~isfield(obj, 'hasHessian'), obj.hasHessian = 0; end
    if ~isfield(obj, 'hasConjHessian'), obj.hasConjHessian = 0; end
end

function prob = ProcessProblem(prob)
    % assign problem attributes to the ones required by minfbe and miname.
    % split the f terms into quadratic and non-quadratic.
    % since we only consider one term for now, assign f to f1 or f2.
    f = ProcessFunction(prob.f);
    prob.g = ProcessFunction(prob.g);
    switch prob.id
        case 1
            if f.isQuadratic
                prob.f1 = f;
                if isfield(prob, 'C'), prob.C1 = prob.C; end
                if isfield(prob, 'Ct'), prob.C1t = prob.Ct; end
                if isfield(prob, 'd'), prob.d1 = prob.d; end
            else
                prob.f2 = f;
                if isfield(prob, 'C'), prob.C2 = prob.C; end
                if isfield(prob, 'Ct'), prob.C2t = prob.Ct; end
                if isfield(prob, 'd'), prob.d2 = prob.d; end
            end
        case 2
            if any(isfield(prob, {'C', 'Ct', 'd'}))
                error('you cannot provide both equality constraints and affine mappings composed with f');
            end
            if f.isConjQuadratic
                prob.f1 = f;
                if isfield(prob, 'A'), prob.A1 = prob.A; end
                if isfield(prob, 'At'), prob.A1t = prob.At; end
            else
                prob.f2 = f;
                if isfield(prob, 'A'), prob.A2 = prob.A; end
                if isfield(prob, 'At'), prob.A2t = prob.At; end
            end
    end
end

function out = ProcessOutput(out)
    % not much to do if id=1.
    % if id=2 must take x1 and x2 and merge them.
    % so far we assume that either f1 or f2 is present (see ProcessProblem)
    % so the output of miname only has either x1 or x2, not both.
    if out.prob.id == 2
        if isfield(out, 'x1')
            out.x = out.x1;
            out = rmfield(out, 'x1');
        elseif isfield(out, 'x2')
            out.x = out.x2;
            out = rmfield(out, 'x2');
        end
    end
end
