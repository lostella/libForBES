function out = miname(prob, opt)
    t0 = tic();
    [prob, dualprob] = ProcessSeparableProblem(prob);
    preprocess = toc(t0);
    if nargin > 1
        dualout = minfbe(dualprob, opt);
    else
        dualout = minfbe(dualprob);
    end
    out = GetPrimalOutput(prob, dualout);
    out.preprocess = preprocess + dualout.preprocess;
    if nargin > 1
        out.opt = opt;
    end
    out.dual = dualout;
end
