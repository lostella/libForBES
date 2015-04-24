function prob = ProcessProblem(prob)
    if ~isfield(prob, 'x0'), error('the starting point x0 must be specified'); end
    if ~isfield(prob, 'useHessian'), prob.useHessian = 0; end
    if ~isfield(prob, 'muf'), prob.muf = 0; end
    prob.n = length(prob.x0);
    prob.Lf = 0;
    eigsOpt.issym = 1;
    eigsOpt.tol = 1e-3;
    if isfield(prob, 'Q') || isfield(prob, 'q')
        prob.istheref1 = true;
        prob.isthereA = true;
        prob.isAfun = false;
        prob.isQfun = false;
        if isfield(prob, 'Q') && isa(prob.Q, 'function_handle')
            prob.isQfun = true;
        elseif ~isfield(prob, 'Q')
            prob.Q = 1;
        end
        if isfield(prob, 'A')
            if isa(prob.A, 'function_handle')
                prob.m1 = length(prob.A(prob.x0));
                if ~isfield(prob, 'AT') || ~isa(prob.AT, 'function_handle')
                    error('must specify both A and AT as function handles');
                end
                prob.isAfun = true;
                if prob.isQfun, funHessian = @(x) prob.AT(prob.Q(prob.A(x)));
                else funHessian = @(x) prob.AT(prob.Q*prob.A(x)); end
            else
                prob.m1 = size(prob.A, 1);
                if prob.isQfun, funHessian = @(x) prob.A'*(prob.Q(prob.A*x));
                else funHessian = @(x) prob.A'*(prob.Q*(prob.A*x)); end
            end
        elseif ~isfield(prob, 'A')
            prob.m1 = prob.n;
            prob.isthereA = false;
            if prob.isQfun, funHessian = @(x) prob.Q(x);
            else funHessian = @(x) prob.Q*x; end
        end
        if isfield(prob, 'Lf1'), prob.Lf = prob.Lf + prob.Lf1;
        else prob.Lf = prob.Lf + eigs(funHessian, prob.n, 1, 'LM', eigsOpt); end
        prob.unknownLf = 0;
        if ~isfield(prob, 'b'), prob.b = zeros(prob.m1, 1); end
        if ~isfield(prob, 'q'), prob.q = zeros(prob.m1, 1); end
    else
        prob.istheref1 = false;
    end
    if isfield(prob, 'f2')
        prob.istheref2 = true;
        prob.isthereC = true;
        prob.isCfun = false;
        if isfield(prob, 'C')
            if isa(prob.C, 'function_handle')
                prob.m2 = length(prob.C(prob.x0));
                if ~isfield(prob, 'CT') || ~isa(prob.CT, 'function_handle')
                    error('must specify both C and CT as function handles');
                end
                prob.isCfun = true;
                funCTC = @(x) prob.CT(prob.C(x));
            else
                prob.m2 = size(prob.C, 1);
                funCTC = @(x) prob.C'*(prob.C*x);
            end
        elseif ~isfield(prob, 'C')
            prob.m2 = prob.n;
            prob.isthereC = false;
            prob.normC = 1;
        end
        if isfield(prob, 'Lf2') && isfield(prob, 'normC')
            prob.Lf = prob.Lf + prob.Lf2*prob.normC^2;
            prob.unknownLf = 0;
        elseif ~isfield(prob, 'Lf2')
            prob.Lf = prob.Lf + 1e-3;
            prob.unknownLf = 1;
        else
            prob.Lf = prob.Lf + prob.Lf2*eigs(funCTC, prob.n, 1, 'LM', eigsOpt);
            prob.unknownLf = 0;
        end
        if ~isfield(prob, 'd'), prob.d = zeros(prob.m2, 1); end
    else
        prob.istheref2 = false;
    end
    if isfield(prob, 'l')
        prob.istherelin = true;
    else
        prob.istherelin = false;
    end
    if prob.istheref1 == false && prob.istheref2 == false, error('missing info'); end
    if ~isfield(prob, 'g'), error('must specify nonsmooth term g'); end
end
