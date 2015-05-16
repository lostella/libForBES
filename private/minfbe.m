function out = minfbe(prob, opt)
    t0 = tic();
    
    if nargin < 1, error('the PROB structure must be provided as first argument'); end
    prob = ProcessCompositeProblem(prob);

    if nargin < 2, opt = []; end
    [opt, name] = ProcessOptions(prob, opt);
    
    lsopt = ProcessLineSearchOptions(opt);
    
    %% initialize output stuff
    objective = zeros(1, opt.maxit);
    ts = zeros(1, opt.maxit);
    residual = zeros(1, opt.maxit);
    msgTerm = '';
    
    %      Q,C1,C2,f2, g
    cnt = [0, 0, 0, 0, 0];

    %% initialize stuff
    x = prob.x0;
    y = x;
    recache = true;
    gam = SelectGamma(prob, opt);
    
    %% initialize specific stuff for the methods
    if opt.method == 2
        S = zeros(prob.n, opt.memory);
        Y = zeros(prob.n, opt.memory);
        YS = zeros(opt.memory, 1);
        LBFGS_col = 1;
        LBFGS_mem = 0;
    end
    if opt.linesearch == 2 % nonmonotone Armijo
        FBErefs = -inf*ones(lsopt.M, 1);
    end
    Q = 0; C = 0; % parameters for Hager-Zhang line search
    
    %% display stuff
    if opt.display >= 2
        fprintf('%6s%11s%11s%11s%11s%11s%11s\n', 'iter', 'gamma', 'optim.', 'object.', '||dir||', 'slope', 'tau');
    end
    
    rejCount = 0;
    flagChangedGamma = 0;
    
    preprocess = toc(t0);
    t0 = tic();
    
    %% main iteration
    for it = 1:opt.maxit
        %% compute proximal gradient step and FBE (if needed)
        if recache || mod(it, opt.recache) == 0
            [cache_current, cnt1] = CacheFBE(prob, gam, y);
            cnt = cnt+cnt1;
            z = cache_current.z;
            if opt.fast && it > 1 && opt.monotone > 0 && ~flagChangedGamma
                if cache_current.FBE > cache_tau.FBE
                    y = cache_tau.z;
                    [cache_current, cnt1] = CacheFBE(prob, gam, y);
                    cnt = cnt+cnt1;
                    rejCount = rejCount+1;
                end
            end
        end
        
        %% check if the Lipschitz constant is to be adjusted
        if opt.adaptive || prob.unknownLf
            flagChangedGamma = 0;
            [fz, cnt1] = Evaluatef(prob, cache_current.z);
            cnt = cnt+cnt1;
            % increase Lf until a candidate Lipschitz constant is found
            while fz + cache_current.gz > cache_current.FBE
                prob.Lf = prob.Lf*2;
                gam = SelectGamma(prob, opt);
                flagChangedGamma = 1;
                [cache_current, cnt1] = CacheFBE(prob, gam, cache_current.x, cache_current);
                [fz, cnt2] = Evaluatef(prob, cache_current.z);
                cnt = cnt+cnt1+cnt2;
            end
        end
        
        %% trace stuff
        ts(1, it) = toc(t0);
        residual(1, it) = norm(cache_current.diff, inf)/gam;
        objective(1, it) = cache_current.FBE;
        
        %% check for termination
        if ~opt.customTerm
            % From sec. 8.2.3.2 of Gill, Murray, Wright (1982).
            absFBE = abs(cache_current.FBE);
            if residual(1, it) <= 10*sqrt(eps) || ...
                    (it > 1 && residual(1, it) <= nthroot(opt.tolOpt, 3)*(1+absFBE) && ...
                    norm(cache_previous.x-cache_current.x, inf) < sqrt(opt.tolOpt)*(1+norm(cache_current.x, inf)) && ...
                    cache_previous.FBE-cache_current.FBE < opt.tolOpt*(1+absFBE))
                msgTerm = [msgTerm, 'reached optimum up to tolOpt'];
                flagTerm = 0;
                break;
            end
        else
            if opt.term(cache_current)
                msgTerm = [msgTerm, 'reached optimum (custom criterion)'];
                flagTerm = 0;
                break;
            end
        end
        
        %% compute gradient of the FBE
        [cache_current, cnt1] = CacheGradFBE(prob, gam, cache_current);
        cnt = cnt+cnt1;
        
        %% compute search direction and slope
        switch opt.method
            case {1, 6}
                dir = -cache_current.gradFBE;
            case 2
                if it == 1 || flagChangedGamma
                    dir = -cache_current.gradFBE; % use steepest descent direction initially
                    LBFGS_col = 1;
                    LBFGS_mem = 0;
                else
                    Sk = cache_current.x - cache_previous.x;
                    Yk = cache_current.gradFBE - cache_previous.gradFBE;
                    YSk = Yk'*Sk;
                    if YSk > 0
                        LBFGS_col = 1+mod(LBFGS_col, opt.memory);
                        LBFGS_mem = min(LBFGS_mem+1, opt.memory);
                        S(:,LBFGS_col) = Sk;
                        Y(:,LBFGS_col) = Yk;
                        YS(LBFGS_col) = YSk;
                    end
                    if LBFGS_mem > 0
                        H = YS(LBFGS_col)/(Y(:,LBFGS_col)'*Y(:,LBFGS_col));
                        dir = LBFGS(S, Y, YS, H, -cache_current.gradFBE, int32(LBFGS_col), int32(LBFGS_mem));
                    else
                        dir = -cache_current.gradFBE;
                    end
                end
            case 3 % CG-DESCENT
                if it == 1 || flagChangedGamma
                    dir = -cache_current.gradFBE; % Initially use steepest descent direction
                else
                    yy = cache_current.gradFBE-cache_previous.gradFBE;
                    dy = dir'*yy;
                    lambda = 1; %Hager-Zhang proposed lambda = 2 but Kou, Dai found that lambda = 1 is more efficient
                    %                 lambda = 2-(dir'*yy)/((dir'*dir)*(yy'*yy));
                    beta = ((yy-lambda*dir*(yy'*yy)/dy)'*cache_current.gradFBE)/dy;
                    etak = -1/(norm(dir)*min(0.01,norm(cache_current.gradFBE)));
                    beta = max(beta,etak);
                    dir = -cache_current.gradFBE + beta*dir;
                end
            case 4 % PRP-CG
                if it == 1 || flagChangedGamma
                    dir = -cache_current.gradFBE; % Initially use steepest descent direction
                else
                    yy = cache_current.gradFBE - cache_previous.gradFBE;
                    beta = max((cache_current.gradFBE'*yy)/(cache_previous.gradFBE'*cache_previous.gradFBE),0);
                    dir = -cache_current.gradFBE + beta*dir;
                end
            case 5 % DYHS-CG
                if it == 1 || flagChangedGamma
                    dir = -cache_current.gradFBE; % Initially use steepest descent direction
                else
                    yy = cache_current.gradFBE - cache_previous.gradFBE;
                    betaDY = (cache_current.gradFBE'*cache_current.gradFBE)/(dir'*yy);
                    betaHS = (cache_current.gradFBE'*yy)/(dir'*yy);
                    beta = max(0,min(betaHS,betaDY));
                    dir = -cache_current.gradFBE + beta*dir;
                end
        end
        slope = cache_current.gradFBE'*dir;
        
        %% precompute stuff for the line search
        [cache_current, cnt1] = CacheLSData(prob, dir, cache_current);
        cnt = cnt+cnt1;

        %% set initial guess for the step length
        switch opt.method
            case 2
                lsopt.tau0 = 1.0;
            case 6
                if it == 1 || flagChangedGamma
                    tau0 = 1.0/norm(cache_current.gradFBE, inf);
                else
                    Sk = cache_current.x-cache_previous.x;
                    Yk = cache_current.gradFBE-cache_previous.gradFBE;
                    tau0 = (Sk'*Sk)/(Sk'*Yk);
                end
                lsopt.tau0 = tau0;
            otherwise
                if it == 1 || flagChangedGamma
                    xinf = norm(cache_current.x,inf);
                    if xinf ~= 0
                        lsopt.tau0 = lsopt.psi0*xinf/norm(cache_current.gradFBE, inf); % g is the gradient at x
                    elseif cache_current.FBE ~= 0
                        % Check changes in Version 5.2 (3). psi0 is set equal
                        % to 0 when x=0;
                        lsopt.tau0 = lsopt.psi0*abs(cache_current.FBE)/(cache_current.gradFBE'*cache_current.gradFBE);
                    else
                        lsopt.tau0 = 1;
                    end
                else
                    %                             lsopt.tau0 = lsopt.psi2*tau;
                    %                             lsopt.tau0 = 2*(cache.FBE - FBE_old)/slope;
                    lsopt.tau0 = -2*max(cache_previous.FBE-cache_current.FBE ,10*eps)/slope;% Fletcher, pp. 38
                    if lsopt.quadStep
                        tp = tau*lsopt.psi1;
                        [cache_tau, cnt1] = DirFBE(prob, gam, tp, cache_current, 1);
                        cnt = cnt+cnt1;
                        if cache_tau.FBE <= cache_current.FBE
                            % quadratic interpolation
                            q = cache_tau.FBE - cache_current.FBE - tp*slope;
                            if q > 0
                                lsopt.tau0 = -(slope*tp^2)/(2*q);
                            end
                        end
                    end
                end
        end

        %% perform the line search
        switch opt.linesearch
            case 1 % backtracking until Armijo condition is met
                [cache_tau, tau, cntLS, info] = ArmijoLS(prob, gam, cache_current, slope, lsopt);
            case 2 % Nonmonotone Armijo backtracking
                FBErefs(1+mod(it, lsopt.M), 1) = cache_current.FBE;
                [cache_tau, tau, cntLS, info] = ArmijoLS(prob, gam, cache_current, slope, lsopt, max(FBErefs));
            case 3 % Lemarechal line search
                [cache_tau, tau, cntLS, info] = LemarechalLS(prob, gam, cache_current, slope, lsopt);
            case 4 % Hager-Zhang line search
                Q = 1 + Q*lsopt.Delta;
                C = C + (abs(cache_current.FBE) - C)/Q;
                [cache_tau, tau, cntLS, info] = HagerZhangLS(prob, gam, cache_current, slope, lsopt);   
                if ~opt.global && ~opt.fast && info ~= 0
                    flagTerm = 2;
                    msgTerm = [msgTerm, 'hager-zhang line search failed at it. ', num2str(it)];
                    break;
                end
                if ~lsopt.AWolfe
                    if abs(cache_tau.FBE-cache_current.FBE) <= lsopt.omega*C % switch to approximate Wolfe conditions
                        lsopt.AWolfe = true;
                    end
                end
            case 5 % More-Thuente line search
                [cache_tau, tau, cntLS, info] = MoreThuenteLS(prob, gam, cache_current, slope, lsopt);
            case 6 % Fletcher's line search (modified matlab function)
                [cache_tau, tau, cntLS, info] = FletcherLS(prob, gam, cache_current, slope,lsopt);
        end
        cnt = cnt+cntLS;

        %% check for line search fails
        if ~opt.global && ~opt.fast && info ~= 0 && opt.linesearch ~= 4
            msgTerm = ['line search failed at it. ', num2str(it), '; '];
            opt.linesearch = 4;
            tau0 = lsopt.tau0;
            lsopt = ProcessLineSearchOptions(opt);
            lsopt.tau0 = tau0;
            [cache_tau, tau, cntLS, info] = HagerZhangLS(prob, gam, cache_current, slope, lsopt);
            cnt = cnt+cntLS;
            if info ~= 0
                flagTerm = 2;
                msgTerm = [msgTerm, 'hager-zhang line search failed at it. ', num2str(it)];
                break;
            end
            if ~lsopt.AWolfe
                if abs(cache_tau.FBE-cache_current.FBE) <= lsopt.omega*C % switch to approximate Wolfe conditions
                    lsopt.AWolfe = true;
                end
            end
        end
        
        %% update the iterate
        if opt.method == 0
            y = cache_current.z;
            cache_previous = cache_current;
            recache = true;
        elseif opt.global
            y = cache_tau.z;
            cache_previous = cache_current;
            recache = true;
        elseif opt.fast
            theta1 = 2/(it+1);
            theta = 2/(it+2);
            v = x + (1/theta1)*(z-x);
            x = cache_tau.z;
            y = theta*v + (1-theta)*x;
            cache_previous = cache_current;
            recache = true;
        else
            y = cache_tau.x;
            cache_previous = cache_current;
            cache_current = cache_tau;
            recache = false;
        end
        
        %% display stuff
        if opt.display == 1
            if mod(it, 100) == 0
                fprintf('.');
            end
            if mod(it, 4000) == 0
                fprintf('\n');
            end
        elseif opt.display >= 2
            fprintf('%6d %7.4e %7.4e %7.4e %7.4e %7.4e %7.4e %d\n', it, gam, residual(1,it), cache_current.FBE, norm(dir), slope, tau, info);
        end
        
    end
    
    if it == opt.maxit
        flagTerm = 1;
        msgTerm = [msgTerm, 'exceeded maximum iterations'];
    end
    
    %% pack up results
    out.name = name;
    out.message = msgTerm;
    out.flag = flagTerm;
    out.gam = gam;
    out.x = cache_current.z;
    if prob.istheref1, out.res1 = cache_current.res1x; end
    if prob.istheref2, out.res2 = cache_current.res2x; end
    out.iterations = it;
    out.operations.cnt_Q = cnt(1);
    out.operations.cnt_C1 = cnt(2);
    out.operations.cnt_C2 = cnt(3);
    out.operations.cnt_f2 = cnt(4);
    out.operations.cnt_g = cnt(5);
    out.operations.rej_extr = rejCount;
    out.residual = residual(1, 1:it);
    out.objective = objective(1, 1:it);
    out.ts = ts(1, 1:it);
    out.prob = prob;
    out.preprocess = preprocess;
    out.opt = opt;
end

function gam = SelectGamma(prob, opt)
    if opt.method == 0 || opt.fast || opt.global
        gam = 1/prob.Lf;
    else
        gam = 0.95/prob.Lf;
    end
end

function [v, cnt] = Evaluatef(prob, x)
    %      Q,C1,C2,f2, g
    cnt = [0, 0, 0, 0, 0];
    f1x = 0; f2x = 0;
    if prob.istheref1
        if prob.isthereC1
            if prob.isC1fun, C1x = prob.C1(x);
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
function [cache, cnt] = CacheFBE(prob, gam, x, cache)
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
    cache.FBE = cache.fx + cache.gz + ...
        cache.gradfx'*cache.diff + ...
        (0.5/gam)*sqnormdiff;
end

%% compute the gradient of the FBE and store reusable quantities
function [cache, cnt] = CacheGradFBE(prob, gam, cache)
    %      Q,C1,C2,f2, g
    cnt = [0, 0, 0, 0, 0];
    Hdiff = 0;
    if prob.istheref1
        if prob.isthereC1
            if prob.isC1fun, C1diff = prob.C1(cache.diff);
            else C1diff = prob.C1*cache.diff; end
            if prob.isQfun, QC1diff = prob.Q(C1diff);
            else QC1diff = prob.Q*C1diff; end
            if prob.isC1fun, C1tQC1diff = prob.C1t(QC1diff);
            else C1tQC1diff = prob.C1'*QC1diff; end
            cnt(2) = cnt(2)+2;
        else
            if prob.isQfun, C1tQC1diff = prob.Q(cache.diff);
            else C1tQC1diff = prob.Q*cache.diff; end
        end
        cnt(1) = cnt(1)+1;
        Hdiff = Hdiff + C1tQC1diff;
    end
    if prob.istheref2
        if prob.isthereC2
            if prob.isC2fun, C2diff = prob.C2(cache.diff);
            else C2diff = prob.C2*cache.diff; end
            cnt(3) = cnt(3)+1;
        else
            C2diff = cache.diff;
        end
        if prob.useHessian
            HC2diff = cache.Hessf2res2x(C2diff);
        else
            res2xepsdiff = cache.res2x + 1e-100i*C2diff;
            [~, gradf2res2xepsd] = prob.callf2(res2xepsdiff);
            cnt(4) = cnt(4)+1;
            HC2diff = imag(gradf2res2xepsd)/1e-100;
        end
        if prob.isthereC2
            if prob.isC2fun, Hdiff = Hdiff + prob.C2t(HC2diff);
            else Hdiff = Hdiff + (prob.C2'*HC2diff); end
            cnt(3) = cnt(3)+1;
        else
            Hdiff = Hdiff + HC2diff;
        end
    end
    cache.gradFBE = (Hdiff - cache.diff/gam);
end

%% precompute quantities needed for the line search procedures
function [cache, cnt] = CacheLSData(prob, dir, cache)
    %      Q,C1,C2,f2, g
    cnt = [0, 0, 0, 0, 0];
    cache.dir = dir;
    if prob.istheref1
        if prob.isthereC1
            if prob.isC1fun, cache.C1dir = prob.C1(dir);
            else cache.C1dir = prob.C1*dir; end
            if prob.isQfun, cache.QC1dir = prob.Q(cache.C1dir);
            else cache.QC1dir = prob.Q*cache.C1dir; end
            if prob.isC1fun, cache.C1tQC1dir = prob.C1t(cache.QC1dir);
            else cache.C1tQC1dir = prob.C1'*cache.QC1dir; end
            cnt(2) = cnt(2)+2;
        else
            cache.C1dir = dir;
            if prob.isQfun, cache.QC1dir = prob.Q(cache.C1dir);
            else cache.QC1dir = prob.Q*cache.C1dir; end
            cache.C1tQC1dir = cache.QC1dir;
        end
        cnt(1) = cnt(1)+1;
        cache.f1linear = cache.gradf1x'*dir;
        cache.f1quad = cache.C1dir'*cache.QC1dir;
    end
    if prob.istheref2
        if prob.isthereC2
            if prob.isC2fun, cache.C2dir = prob.C2(dir);
            else cache.C2dir = prob.C2*dir; end
            cnt(3) = cnt(3)+1;
        else
            cache.C2dir = dir;
        end
    end
    if prob.istherelin
        cache.lindir = prob.l'*dir;
    end
end
