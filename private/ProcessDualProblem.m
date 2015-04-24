function dualprob = ProcessDualProblem(prob)
    if isfield(prob, 'c')
        if norm(prob.c) > 0, dualprob.l = prob.c; end
        n = length(prob.c);
        dualprob.x0 = zeros(n, 1);
    else
        error('the number right hand side c of the constraint must be specified');
    end
    if not(isfield(prob, 'x1step') || isfield(prob, 'x2step'))
        error('both smooth terms f1 and f2 are missing');
    end
    if isfield(prob, 'x1step')
        if isfield(prob, 'A1')
            if isa(prob.A1, 'function_handle')
                if ~isfield(prob, 'A1T') || ~isa(prob.A1T, 'function_handle')
                    error('must specify both A1 and A1T as function handles');
                end
                prob.isA1fun = true;
                gf1c0 = prob.x1step(-prob.A1T(zeros(n, 1)));
                dualprob.Q = @(w) prob.x1step(w)-gf1c0;
                dualprob.A = @(x) -prob.A1T(x);
                dualprob.AT = @(y) -prob.A1(y);
                dualprob.q = gf1c0;
            else
                prob.isA1fun = false;
                gf1c0 = prob.x1step(-prob.A1'*zeros(n, 1));
                dualprob.Q = @(w) prob.x1step(w)-gf1c0;
                dualprob.A = -prob.A1';
                dualprob.q = gf1c0;
            end
        else
            error('must specify matrix A1 in the constraint');
        end
        if isfield(prob, 'muf1'), dualprob.Lf1 = 1/prob.muf1; end
    end
    if isfield(prob, 'x2step')
        dualprob.f2 = @(y) ConjugateFunction(prob.x2step, y);
        if isfield(prob, 'A2')
            if isa(prob.A2, 'function_handle')
                if ~isfield(prob, 'A2T') || ~isa(prob.A2T, 'function_handle')
                    error('must specify both A2 and A2T as function handles');
                end
                prob.isA2fun = true;
                dualprob.C = -prob.A2T;
                dualprob.CT = -prob.A2;
            else
                prob.isA2fun = false;
                dualprob.C = -prob.A2';
            end
        else
            error('must specify matrix A2 in the constraint');
        end
        if isfield(prob, 'muf2'), dualprob.Lf2 = 1/prob.muf2; end
        if isfield(prob, 'normA2'), dualprob.normC = prob.normA2; end
    end
    if isfield(prob, 'B')
        if isa(prob.B, 'function_handle'), error('must specify matrix B *as a matrix* in the constraint'); end
    else
        error('must specify matrix B in the constraint');
    end
    dualprob.g = @(y, gam) ProxConjugateFunction(prob.zstep, prob.B, y, gam);
end

function [val, x] = ConjugateFunction(xstep, w)
% Argument xstep computes x and v, given y
%
%   x = argmin{f(w) - <y, w>}
%   v = f(x), where x is defined as above.
%
% Therefore x and w are conjugate points with respect to f and f*,
% and by the conjugate subgradient theorem one has
%
%   f*(w) = <w,x> - f(x).
%
    [x, fx] = xstep(y);
    val = w'*x-fx;
end

function [prox, val] = ProxConjugateFunction(zstep, B, y, gam)
% Argument zstep(y, gam) computes z and v:
%   
%   z = argmin{g(z) + <y, Bz> + gam/2 ||Bz||^2}
%   v = g(z), where z is defined as above.
%
% Equivalently
%
%   z = argmin{g(z) + gam/2 ||y/gam + Bz||^2}
%
% By the optimality conditions for z one has
%
%   -B'(y + gam Bz) is a subgradient of g at z,
%   z is a subgradient of g* at -B'(y + gam Bz)
%   (by conjugate subgradient theorem) and
%   g*(-B'(y + gam Bz)) = -z'B'(y + gam Bz) - g(z)
%
    [z, v] = zstep(y, gam);
    Bz = B*z;
    prox = y+gam*Bz;
    val = -prox'*Bz - v;
end
