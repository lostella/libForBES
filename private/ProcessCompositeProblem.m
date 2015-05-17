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

function prob = ProcessCompositeProblem(prob)
    if ~isfield(prob, 'x0'), error('the starting point x0 must be specified'); end
    if ~isfield(prob, 'useHessian'), prob.useHessian = 0; end
    if ~isfield(prob, 'muf'), prob.muf = 0; end
    prob.n = length(prob.x0);
    prob.Lf = 0;
    eigsOpt.issym = 1;
    eigsOpt.tol = 1e-3;
    if any(isfield(prob, {'Q', 'q'}))
        prob.istheref1 = true;
        prob.isthereC1 = true;
        prob.isC1fun = false;
        prob.isQfun = false;
        if isfield(prob, 'Q') && isa(prob.Q, 'function_handle')
            prob.isQfun = true;
        elseif ~isfield(prob, 'Q')
            prob.Q = 1;
        end
        if isfield(prob, 'C1')
            if isa(prob.C1, 'function_handle')
                prob.m1 = length(prob.C1(prob.x0));
                if ~isfield(prob, 'C1t') || ~isa(prob.C1t, 'function_handle')
                    error('must specify both C1 and C1t as function handles');
                end
                prob.isC1fun = true;
                if prob.isQfun, funHessian = @(x) prob.C1t(prob.Q(prob.C1(x)));
                else funHessian = @(x) prob.C1t(prob.Q*prob.C1(x)); end
            else
                prob.m1 = size(prob.C1, 1);
                if prob.isQfun, funHessian = @(x) prob.C1'*(prob.Q(prob.C1*x));
                else funHessian = @(x) prob.C1'*(prob.Q*(prob.C1*x)); end
            end
        else
            prob.m1 = prob.n;
            prob.isthereC1 = false;
            if prob.isQfun, funHessian = @(x) prob.Q(x);
            else funHessian = @(x) prob.Q*x; end
        end
        if isfield(prob, 'Lf1'), prob.Lf = prob.Lf + prob.Lf1;
        else prob.Lf = prob.Lf + eigs(funHessian, prob.n, 1, 'LM', eigsOpt); end
        prob.unknownLf = 0;
        if ~isfield(prob, 'd1'), prob.d1 = zeros(prob.m1, 1); end
        if ~isfield(prob, 'q'), prob.q = zeros(prob.m1, 1); end
    else
        prob.istheref1 = false;
    end
    if isfield(prob, 'f2')
        if ~isfield(prob.f2, 'makef'), error('function of f2 is not defined'); end
        prob.callf2 = prob.f2.makef();
        prob.istheref2 = true;
        prob.isthereC2 = true;
        prob.isC2fun = false;
        if isfield(prob, 'C2')
            if isa(prob.C2, 'function_handle')
                prob.m2 = length(prob.C2(prob.x0));
                if ~isfield(prob, 'C2t') || ~isa(prob.C2t, 'function_handle')
                    error('must specify both C2 and C2t as function handles');
                end
                prob.isC2fun = true;
                funC2tC2 = @(x) prob.C2t(prob.C2(x));
            else
                prob.m2 = size(prob.C2, 1);
                funC2tC2 = @(x) prob.C2'*(prob.C2*x);
            end
        else
            prob.m2 = prob.n;
            prob.isthereC2 = false;
            prob.normC2 = 1;
        end
        if isfield(prob.f2, 'L') && isfield(prob, 'normC2')
            prob.Lf = prob.Lf + prob.f2.L*prob.normC2^2;
            prob.unknownLf = 0;
        elseif ~isfield(prob.f2, 'L')
            x1 = prob.x0+1e-6;
            if prob.isC2fun, C2x0 = prob.C2(prob.x0); else C2x0 = prob.C2*prob.x0; end
            if prob.isC2fun, C2x1 = prob.C2(x1); else C2x1 = prob.C2*x1; end
            [~, gradf0] = prob.callf2(C2x0);
            [~, gradf1] = prob.callf2(C2x1);
            prob.Lf = prob.Lf + norm(gradf0-gradf1)/1e-6;
            prob.unknownLf = 1;
        else
            prob.Lf = prob.Lf + prob.f2.L*eigs(funC2tC2, prob.n, 1, 'LM', eigsOpt);
            prob.unknownLf = 0;
        end
        if ~isfield(prob, 'd2'), prob.d2 = zeros(prob.m2, 1); end
    else
        prob.istheref2 = false;
    end
    if isfield(prob, 'l')
        prob.istherelin = true;
    else
        prob.istherelin = false;
    end
    if prob.istheref1 == false && prob.istheref2 == false, error('you must specify at least one of f1 and f2'); end
    if ~isfield(prob, 'g'), error('you must specify the nonsmooth term g'); end
    if ~isfield(prob.g, 'makeprox'), error('the prox for the term g you specified is not available'); end
    prob.callg = prob.g.makeprox();
end
