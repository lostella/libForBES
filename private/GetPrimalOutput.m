function out = GetPrimalOutput(prob, dualout)
    out.name = dualout.name;
    out.message = dualout.message;
    out.flag = dualout.flag;
    out.gam = dualout.gam;
    Ax = 0;
    out.y = dualout.x;
    if isfield(prob, 'x1step')
        if isa(prob.A1, 'function_handle')
            out.x1 = prob.x1step(-prob.A1T(out.y));
            Ax = Ax+prob.A1(out.x1);
        else
            out.x1 = prob.x1step(-prob.A1'*out.y);
            Ax = Ax+prob.A1*out.x1;
        end
    end
    if isfield(prob, 'x2step')
        if isa(prob.A2, 'function_handle')
            out.x2 = prob.x2step(-prob.A2T(out.y));
            Ax = Ax+prob.A2(out.x2);
        else
            out.x2 = prob.x2step(-prob.A2'*out.y);
            Ax = Ax+prob.A2*out.x2;
        end
    end
    if isfield(prob, 'c')
        [out.z, ~] = prob.zstep(out.y+out.gam*(Ax-prob.c), out.gam);
    else
        [out.z, ~] = prob.zstep(out.y+out.gam*Ax, out.gam);
    end
    out.iterations = dualout.iterations;
    if isfield(dualout, 'operations'), out.operations = dualout.operations; end
    out.residual = dualout.residual;
    out.ts = dualout.ts;
    out.preprocess = dualout.preprocess;
    out.prob = prob;
    out.dual = dualout;
end
