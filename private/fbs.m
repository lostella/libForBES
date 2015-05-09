function out = fbs(prob, opt)

    t0 = tic();
    
    if nargin < 1, error('the PROB structure must be provided as first argument'); end
    prob = ProcessCompositeProblem(prob);
    
    if ~isfield(opt, 'tolOpt'), opt.tolOpt = 1e-5; end
    if ~isfield(opt, 'tolProg'), opt.tolProg = 0; end
    if ~isfield(opt, 'term'), opt.customTerm = false;
    else opt.customTerm = true; end
    if ~isfield(opt, 'maxit'), opt.maxit = 10*prob.n; end
    if ~isfield(opt, 'fast'), opt.fast = 1; end
    if ~isfield(opt, 'monotone'), opt.monotone = 0; end
    if ~isfield(opt, 'display'), opt.display = 1; end
    
    %% initialize output stuff
    residual = zeros(1, opt.maxit);
    ts = zeros(1, opt.maxit);
    objective = zeros(1, opt.maxit);
    msgTerm = '';
    
    gam = 1/prob.Lf;
    xk = prob.x0;
    vk = prob.x0;
    
    preprocess = toc(t0);
    t0 = tic();
    
    for it = 1:opt.maxit
        if opt.fast
            if prob.muf == 0
                theta = 2/(it+1); % since it starts from 1
            else
                theta = sqrt(prob.muf/prob.Lf);
            end
            yk = (1-theta)*xk+theta*vk;
        else
            yk = xk;
        end
        
        cache_yk = ForwardBackwardStep(prob, gam, yk);
        u = cache_yk.z;
        
        ts(1, it) = toc(t0);
        residual(1, it) = norm(cache_yk.diff, inf)/gam;
        objective(1, it) = cache_yk.fx + cache_yk.gz + ...
            cache_yk.gradfx'*cache_yk.diff + ...
            (0.5/gam)*cache_yk.normdiff^2;
        
        if opt.customTerm
            if opt.term(cache_yk)
                msgTerm = [msgTerm, 'reached optimum (custom criterion)'];
                flagTerm = 0;
                break;
            end
        else
            if residual(1, it) <= opt.tolOpt
                msgTerm = [msgTerm, 'reached optimum up to tolOpt'];
                flagTerm = 0;
                break;
            end
        end
        
        xknew = u;
        
        if opt.fast
            vk = xk + (cache_yk.z-xk)/theta;
        end
        xk = xknew;
        
        if opt.display == 1
            if mod(it, 100) == 0
                fprintf('.');
            end
            if mod(it, 4000) == 0
                fprintf('\n');
            end
        end
    end
    
    if it == opt.maxit
        flagTerm = 1;
        msgTerm = [msgTerm, 'exceeded maximum iterations'];
    end
    
    out.name = ['FBS, fast=', num2str(opt.fast), ', mon=', num2str(opt.monotone)'];
    out.message = msgTerm;
    out.flag = flagTerm;
    out.gam = gam;
    out.x = cache_yk.x;
    out.z = cache_yk.z;
    if prob.istheref1, out.res1 = cache_yk.res1x; end
    if prob.istheref2, out.res2 = cache_yk.res2x; end
    out.iterations = it;
    out.operations.cnt_A = 2*it;
    out.operations.cnt_C = 2*it;
    out.operations.cnt_f2 = it;
    out.residual = residual(1, 1:it);
    out.objective = objective(1, 1:it);
    out.ts = ts(1, 1:it);
    out.prob = prob;
    out.preprocess = preprocess;
    out.opt = opt;
end

function [v, cnt] = Evalf(prob, x)
    %      Q,C1,C2,f2, g
    cnt = [0, 0, 0, 0, 0];
    f1x = 0; f2x = 0;
    if prob.istheref1
        if prob.isthereC1
            if prob.isAfun, C1x = prob.C1(x);
            else C1x = prob.C1*x; end
            res1x = C1x - prob.d1;
            if prob.isQfun, Qres1x = prob.Q(res1x);
            else Qres1x = prob.Q*res1x; end
            cnt(2) = cnt(2)+1;
        else
            res1x = x - prob.d1;
            if prob.isQfun, Qres1x = prob.Q(res1x);
            else Qres1x = prob.Q*res1x; end
        end
        cnt(1) = cnt(1)+1;
        f1x = 0.5*(res1x'*Qres1x) + prob.q'*res1x;
    end
    if prob.istheref2
        if prob.isthereC2
            if prob.isC2fun, C2x = prob.C2(x);
            else C2x = prob.C2*x; end
            res2x = C2x - prob.d2;
            f2x = prob.callf2(res2x);
            cnt(3) = cnt(3)+1;
        else
            res2x = x - prob.d2;
            f2x = prob.callf2(res2x);
        end
        cnt(4) = cnt(4)+1;
    end
    if prob.istherelin
        v = f1x + f2x + prob.l'*x;
    else
        v = f1x + f2x;
    end
end

%% compute the FBE value and store reusable quantities
function [cache, cnt] = ForwardBackwardStep(prob, gam, x, cache)
    %      Q,C1,C2,f2, g
    cnt = [0, 0, 0, 0, 0];
    f1x = 0; gradf1x = 0;
    f2x = 0; gradf2x = 0;
    if nargin < 4
        cache.x = x;
        if prob.istheref1
            if prob.isthereC1
                if prob.isC1fun, C1x = prob.C1(cache.x);
                else C1x = prob.C1*cache.x; end
                cache.res1x = C1x - prob.d1;
                if prob.isQfun, cache.Qres1x = prob.Q(cache.res1x);
                else cache.Qres1x = prob.Q*cache.res1x; end
                if prob.isC1fun, gradf1x = prob.C1t(cache.Qres1x + prob.q);
                else gradf1x = prob.C1'*(cache.Qres1x + prob.q); end
                cnt(2) = cnt(2)+2;
            else
                cache.res1x = cache.x - prob.d1;
                if prob.isQfun, cache.Qres1x = prob.Q(cache.res1x);
                else cache.Qres1x = prob.Q*cache.res1x; end
                gradf1x = cache.Qres1x + prob.q;
            end
            cnt(1) = cnt(1)+1;
            cache.gradf1x = gradf1x;
            f1x = 0.5*(cache.res1x'*cache.Qres1x) + prob.q'*cache.res1x;
            cache.f1x = f1x;
        end
        if prob.istheref2
            if prob.isthereC2
                if prob.isC2fun, C2x = prob.C2(cache.x);
                else C2x = prob.C2*cache.x; end
                cache.res2x = C2x - prob.d2;
                if prob.useHessian
                    [f2x, gradf2res2x, cache.Hessf2res2x] = prob.callf2(cache.res2x);
                else
                    [f2x, gradf2res2x] = prob.callf2(cache.res2x);
                end
                if prob.isC2fun, gradf2x = prob.C2t(gradf2res2x);
                else gradf2x = prob.C2'*gradf2res2x; end
                cnt(3) = cnt(3)+2;
            else
                cache.res2x = cache.x - prob.d2;
                if prob.useHessian
                    [f2x, gradf2res2x, cache.Hessf2res2x] = prob.callf2(cache.res2x);
                else
                    [f2x, gradf2res2x] = prob.callf2(cache.res2x);
                end
                gradf2x = gradf2res2x;
            end
            cnt(4) = cnt(4)+1;
            cache.gradf2x = gradf2x;
        end
        if prob.istherelin
            cache.flinx = prob.l'*x;
            cache.fx = f1x + f2x + cache.flinx;
            cache.gradfx = gradf1x + gradf2x + prob.l;
        else
            cache.fx = f1x + f2x;
            cache.gradfx = gradf1x + gradf2x;
        end
    end
    y = cache.x - gam*cache.gradfx;
    [cache.z, cache.gz] = prob.callg(y, gam);
    cnt(5) = cnt(5)+1;
    cache.diff = cache.z-cache.x;
    sqnormdiff = cache.diff'*cache.diff;
    cache.normdiff = sqrt(sqnormdiff);
end

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
                if ~isfield(prob, 'C1t') || ~isa(prob.AT, 'function_handle')
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
                prob.isCfun = true;
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
            prob.Lf = prob.Lf + 1e-3;
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
