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
        fprintf('%6s%11s%11s%11s%11s%11s\n', 'iter', 'gamma', 'optim.', '||dir||', 'slope', 'tau');
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
        if opt.customTerm
            if opt.term(cache_current)
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
        
        %% compute gradient of the FBE (not needed by FBS)
        if opt.method > 0
            [cache_current, cnt1] = CacheGradFBE(prob, gam, cache_current);
            cnt = cnt+cnt1;
        end
        
        %% compute search direction and slope
        switch opt.method
            case {1, 6}
                dir = -cache_current.gradFBE;
                slope = cache_current.gradFBE'*dir;
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
                slope = cache_current.gradFBE'*dir;
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
                slope = cache_current.gradFBE'*dir;
            case 4 % PRP-CG
                if it == 1 || flagChangedGamma
                    dir = -cache_current.gradFBE; % Initially use steepest descent direction
                else
                    yy = cache_current.gradFBE - cache_previous.gradFBE;
                    beta = max((cache_current.gradFBE'*yy)/(cache_previous.gradFBE'*cache_previous.gradFBE),0);
                    dir = -cache_current.gradFBE + beta*dir;
                end
                slope = cache_current.gradFBE'*dir;
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
                slope = cache_current.gradFBE'*dir;
        end
        
        %% perform line search
        % precompute other stuff for the line search
        [cache_current, cnt1] = CacheLSData(prob, dir, cache_current);
        cnt = cnt+cnt1;

        % set initial guess for the step length
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
                % test wether we enforce the Wolfe or approximate Wolfe
                if lsopt.PertRule
                    lsopt.epsilon = lsopt.eps*abs(cache_current.FBE);
                else
                    lsopt.epsilon = lsopt.eps;
                end
                [cache_tau, tau, cntLS, info] = HagerZhangLS(prob, gam, cache_current, slope, lsopt);   
                if ~opt.global && info ~= 0
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

        % check for line search fails
        if ~opt.global && ~opt.fast && info ~= 0 && opt.linesearch ~= 4
            msgTerm = ['line search failed at it. ', num2str(it), '; '];
            opt.linesearch = 4;
            tau0 = lsopt.tau0;
            lsopt = ProcessLineSearchOptions(opt);
            lsopt.tau0 = tau0;
            if lsopt.PertRule
                lsopt.epsilon = lsopt.eps*abs(cache_current.FBE);
            else
                lsopt.epsilon = lsopt.eps;
            end
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
            fprintf('%6d %7.4e %7.4e %7.4e %7.4e %7.4e\n', it, gam, residual(1,it), norm(dir), slope, tau);
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

function [opt, name] = ProcessOptions(prob, opt)
    % fill in missing options with defaults
    if ~isfield(opt, 'tolOpt'), opt.tolOpt = 1e-5; end
    if ~isfield(opt, 'term'), opt.customTerm = false;
    else opt.customTerm = true; end
    if ~isfield(opt, 'maxit'), opt.maxit = 10*prob.n; end
    if ~isfield(opt, 'method'), opt.method = 'lbfgs'; end
    if ~isfield(opt, 'linesearch')
        switch opt.method
            case 'sd'
                opt.linesearch = 'armijo';
            case 'lbfgs'
                opt.linesearch = 'hager-zhang';
            case 'cg-desc'
                opt.linesearch = 'hager-zhang';
            case 'cg-prp'
                opt.linesearch = 'hager-zhang';
            case 'cg-dyhs'
                opt.linesearch = 'hager-zhang';
            case 'bb'
                opt.linesearch = 'nonmonotone-armijo';
        end
    end
    if ~isfield(opt, 'variant'), opt.variant = 'global'; end
    if ~isfield(opt, 'recache'), opt.recache = 100; end
    if ~isfield(opt, 'memory'), opt.memory = 11; end
    if ~isfield(opt, 'adaptive'), opt.adaptive = 0; end
    if ~isfield(opt, 'display'), opt.display = 0; end
    name = [opt.method,', ', opt.linesearch, ', ', opt.variant];
    % translate labels into integer codes
    switch opt.method
        case 'sd'
            opt.method = 1;
        case 'lbfgs'
            opt.method = 2;
        case 'cg-desc'
            opt.method = 3;
        case 'cg-prp'
            opt.method = 4;
        case 'cg-dyhs'
            opt.method = 5;
        case 'bb'
            opt.method = 6;
        otherwise
            error('unknown method');
    end
    switch opt.linesearch
        case 'armijo'
            opt.linesearch = 1;
        case 'nonmonotone-armijo'
            opt.linesearch = 2;
        case 'lemarechal'
            opt.linesearch = 3;
        case 'hager-zhang'
            opt.linesearch = 4;
        case 'more-thuente'
            opt.linesearch = 5;
        case 'fletcher'
            opt.linesearch = 6;
        otherwise
            error('unknown line search');
    end
    switch opt.variant
        case 'basic'
            opt.fast = 0;
            opt.global = 0;
            opt.monotone = 0;
        case 'global'
            opt.fast = 0;
            opt.global = 1;
            opt.monotone = 0;
        case 'fast'
            opt.fast = 1;
            opt.global = 0;
            opt.monotone = 1;
        otherwise
            error('unknown variant');
    end
end

function lsopt = ProcessLineSearchOptions(opt)
    %  factor in [0, 1] used to compute average cost magnitude C_k as follows:
    % Q_k = 1 + (Delta)Q_k-1, Q_0 = 0,  C_k = C_k-1 + (|f_k| - C_k-1)/Q_k
    lsopt.Delta = 0.7;% this goes here to include Hager-Zhang line search as a backup
    % Wolfe line search parameter delta, range [0, .5]
    % phi (a) - phi (0) <= delta phi'(0)
    lsopt.delta = 0.1;
    switch opt.linesearch
        case 1 % armijo backtracking
%         lsopt.str = 'Armijo Backtracking';
            lsopt.progTol = 0;
            lsopt.nLS = 50;
        case 2 % Nonmonotone Armijo
%         lsopt.str = 'Nonmonotone Armijo Backtracking';
            lsopt.progTol = 0;
            lsopt.nLS = 50;
            lsopt.M = 5;
        case 3 % Lemarechal line search
            lsopt.sigma = 0.9;
            % maximum number of iterations
            lsopt.nbracket = 100;
            % type of interpolation - 0 [bisection], 1 [quadratick
            % interpolation], 2 [cubic interpolation when possible]
            lsopt.interp = 1;
            if isfield(opt, 'interp'), lsopt.interp = opt.interp; end
            % stop when length of interval is below progTol
            lsopt.progTol = 0;
            %  growth factor in search for initial bracket interval
            lsopt.rho = 5;
            % parameter for safe-guarding (must be in (0,1/2])
            lsopt.theta = 0.49;
        case 4 % HagerZhang
            lsopt.sigma = 0.9;
            %  maximum number of times the bracketing interval grows during expansion
            lsopt.nexpand = 50;
            % maximum number of secant steps
            lsopt.nsecant = 50;
            % maximum number of times the bracketing interval contracts
            lsopt.ncontract = 10;
            % factor by which eps grows when line search fails during contraction
            lsopt.egrow = 10;
            lsopt.QuadOK = true;
            % T => when possible, use a cubic step in the line search
            lsopt.UseCubic = true;
            % true => estimated error in function value is eps*Ck,
            % false => estimated error in function value is eps */
            lsopt.PertRule = true;
            lsopt.eps = 1e-6;
            % |f| < SmallCost*starting cost => skip QuadStep and set PertRule = FALSE*/
            lsopt.SmallCost = 1e-30;
            % T => use approximate Wolfe line search
            % F => use ordinary Wolfe line search, switch to approximate Wolfe when
            % |f_k+1-f_k| < omega*C_k, C_k = average size of cost */
            lsopt.AWolfe = false;
            lsopt.omega = 1e-3;
            % factor by which secant step is amplified during expansion phase where minimizer is bracketed
            lsopt.SecantAmp = 1.05;
            % factor by which rho grows during expansion phase where minimizer is bracketed
            lsopt.RhoGrow = 2.0;
            %  maximum number of times that eps is updated
            lsopt.neps = 5;
            % maximum factor secant step increases stepsize in expansion phase
            lsopt.ExpandSafe = 200;
            % value of the parameter theta in the cg_descent update formula:
            % W. W. Hager and H. Zhang, A survey of nonlinear conjugate gradient
            % methods, Pacific Journal of Optimization, 2 (2006), pp. 35-58.
            lsopt.theta = 0.5;
            %  growth factor in search for initial bracket interval
            lsopt.rho = 5;
            % decay factor for bracket interval width in line search, range (0, 1)
            lsopt.gamma = 0.66;
        case 5 % More Thuente
            lsopt.sigma = 0.9;
            lsopt.progTol = 0;
            lsopt.tmin = 0;
            lsopt.tmax = 1e15;
            lsopt.maxfev = 100;
        case 6 % Fletcher
            lsopt.sigma = 0.9;
            %  maximum number of times the bracketing interval grows during expansion
            lsopt.nbracket = 50;
            % maximum number of section steps
            lsopt.nsection = 50;
            % stop when progress is below progTol
            lsopt.progTol = 0;
            % estimate of minimum value of the function
            lsopt.fmin = -inf;
    end

    % if method is not FBS or L-BFGS then initial stepsize is selected
    % according to Hager-Zhang
    if opt.method ~= 2
        lsopt.quadStep = true;
        % starting guess for line search =
        % psi0 ||x_0||_infty over ||g_0||_infty if x_0 != 0
        % psi0 |f(x_0)|/||g_0||_2               otherwise */
        lsopt.psi0 = 0.01;
        % when the function is approximately quadratic, use gradient at
        % psi1*psi2*previous step for estimating initial stepsize */
        lsopt.psi1 = 1.0 ;
        % when starting a new cg iteration, our initial guess for the line
        % search stepsize is psi2*previous step */
        lsopt.psi2 = 2;
    end
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
