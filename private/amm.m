% Solves problems of the form
%
%   minimize f(x) + g(z)
%   subject to Ax + Bz = 0
%
function out = amm(prob, opt)
    t0 = tic();
    [prob, dualprob] = ProcessSeparableProblem(prob);
    preprocess = toc(t0);
    if nargin > 1
        dualout = fbs(dualprob, opt);
    else
        dualout = fbs(dualprob);
    end
    out = GetPrimalOutput(prob, dualout);
    out.preprocess = preprocess + dualout.preprocess;
    if nargin > 1
        out.opt = opt;
    end
    out.dual = dualout;
end

function [prob, dualprob] = ProcessSeparableProblem(prob)
    if isfield(prob, 'b')
        if norm(prob.b) > 0, dualprob.l = prob.b; end
        n = length(prob.b);
        if isfield(prob, 'y0'), dualprob.x0 = prob.y0;
        else dualprob.x0 = zeros(n, 1); end
    else
        error('you must specify the right hand side b of the equality constraint');
    end
    if ~any(isfield(prob, {'f1', 'f2'}))
        error('both smooth terms f1 and f2 are missing');
    end
    if isfield(prob, 'f1')
        if ~isfield(prob.f1, 'makefconj'), error('conjugate function of f1 is not defined'); end
        prob.callf1conj = prob.f1.makefconj();
        if isfield(prob, 'A1')
            if isa(prob.A1, 'function_handle')
                if ~isfield(prob, 'A1t') || ~isa(prob.A1t, 'function_handle')
                    error('you must specify both A1 and A1t as function handles');
                end
                prob.isA1fun = true;
                [dualprob.Q, dualprob.q] = make_quad_conj(prob.f1, prob.A1t(zeros(n, 1)));
                dualprob.C1 = @(x) -prob.A1t(x);
                dualprob.C1t = @(y) -prob.A1(y);
            else
                prob.isA1fun = false;
                [dualprob.Q, dualprob.q] = make_quad_conj(prob.f1, prob.A1'*zeros(n, 1));
                dualprob.C1 = -prob.A1';
            end
        else
            error('you must specify matrix A1 in the constraint');
        end
        if isfield(prob, 'muf1'), dualprob.Lf1 = 1/prob.muf1; end
    end
    if isfield(prob, 'f2')
        if ~isfield(prob.f2, 'makefconj'), error('conjugate function of f2 is not defined'); end
        prob.callf2conj = prob.f2.makefconj();
        dualprob.f2.makef = @() prob.f2.makefconj();
        if isfield(prob, 'A2')
            if isa(prob.A2, 'function_handle')
                if ~isfield(prob, 'A2t') || ~isa(prob.A2T, 'function_handle')
                    error('you must specify both A2 and A2t as function handles');
                end
                prob.isA2fun = true;
                dualprob.C2 = -prob.A2t;
                dualprob.C2t = -prob.A2;
            else
                prob.isA2fun = false;
                dualprob.C2 = -prob.A2';
            end
        else
            error('yout must specify matrix A2 in the constraint');
        end
        if isfield(prob, 'muf2'), dualprob.Lf2 = 1/prob.muf2; end
        if isfield(prob, 'normA2'), dualprob.normC2 = prob.normA2; end
    end
    if isfield(prob, 'B')
        if isa(prob.B, 'function_handle'), error('you must specify matrix B as a matrix in the constraint'); end
    else
        error('you must specify matrix B in the constraint');
    end
    mus = sum((prob.B).*(prob.B), 1);
    if (max(mus)-min(mus))/max(mus) > 10*eps, error('B''*B must be a multiple of the identity'); end
    prob.muB = mus(1);
    if ~isfield(prob, 'g'), error('you must specify term g'); end
    if ~isfield(prob.g, 'makeprox'), error('the prox for the term g you specified is not available'); end
    prob.callg = prob.g.makeprox();
    dualprob.g.makeprox = @() make_prox_conj(prob.callg, prob.B, prob.muB);
end

function op = make_prox_conj(prox, B, mu)
    op = @(y, gam) call_prox_conj(y, gam, prox, B, mu);
end

function [proxpoint, proxval] = call_prox_conj(y, gam, prox, B, mu)
    mugam = mu*gam;
    [z, v] = prox((-B'*y)/mugam, 1/mugam);
    Bz = B*z;
    proxpoint = y+gam*Bz;
    proxval = -proxpoint'*Bz - v;
end

function [Q, q] = make_quad_conj(obj, zero)
    fconj = obj.makefconj();
    [~, q] = fconj(zero);
    Q = @(x) call_quad_conj(x, fconj, q);
end
 
function y = call_quad_conj(x, fconj, gradfconj0)
    [~, gradfconjx] = fconj(x);
    y = gradfconjx-gradfconj0;
end

function out = GetPrimalOutput(prob, dualout)
    out.name = dualout.name;
    out.message = dualout.message;
    out.flag = dualout.flag;
    out.gam = dualout.gam;
    Ax = 0;
    y = dualout.x;
    if isfield(prob, 'f1')
        if isa(prob.A1, 'function_handle')
            [~, out.x1] = prob.callf1conj(-prob.A1t(y));
            Ax = Ax+prob.A1(out.x1);
        else
            [~, out.x1] = prob.callf1conj(-prob.A1'*y);
            Ax = Ax+prob.A1*out.x1;
        end
    end
    if isfield(prob, 'f2')
        if isa(prob.A2, 'function_handle')
            [~, out.x2] = prob.callf2conj(-prob.A2t(y));
            Ax = Ax+prob.A2(out.x2);
        else
            [~, out.x2] = prob.callf2conj(-prob.A2'*y);
            Ax = Ax+prob.A2*out.x2;
        end
    end
    mugam = prob.muB*out.gam;
    if isfield(prob, 'b')
        [out.z, ~] = prob.callg(-prob.B'*(y+out.gam*(Ax-prob.b))/mugam, 1/mugam);
    else
        [out.z, ~] = prob.callg(-prob.B'*(y+out.gam*Ax)/mugam, 1/mugam);
    end
    out.y = y;
    out.iterations = dualout.iterations;
    if isfield(dualout, 'operations'), out.operations = dualout.operations; end
    out.residual = dualout.residual;
    out.ts = dualout.ts;
    out.prob = prob;
    out.preprocess = dualout.preprocess;
end

