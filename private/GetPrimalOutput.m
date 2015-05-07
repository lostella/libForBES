function out = GetPrimalOutput(prob, dualout)
    out.name = dualout.name;
    out.message = dualout.message;
    out.flag = dualout.flag;
    out.gam = dualout.gam;
    Ax = 0;
    out.y = dualout.x;
    if isfield(prob, 'f1')
        if isa(prob.A1, 'function_handle')
            [~, out.x1] = prob.callf1conj(-prob.A1t(out.y));
            Ax = Ax+prob.A1(out.x1);
        else
            [~, out.x1] = prob.callf1conj(-prob.A1'*out.y);
            Ax = Ax+prob.A1*out.x1;
        end
    end
    if isfield(prob, 'f2')
        if isa(prob.A2, 'function_handle')
            [~, out.x2] = prob.callf2conj(-prob.A2t(out.y));
            Ax = Ax+prob.A2(out.x2);
        else
            [~, out.x2] = prob.callf2conj(-prob.A2'*out.y);
            Ax = Ax+prob.A2*out.x2;
        end
    end
    mugam = prob.muB*out.gam;
    if isfield(prob, 'b')
        [out.z, ~] = prob.callg(-prob.B'*(out.y+out.gam*(Ax-prob.b))/mugam, 1/mugam);
%         [out.z, ~] = prob.callg(out.y/out.gam + (Ax-prob.c), 1/out.gam);
    else
        [out.z, ~] = prob.callg(-prob.B'*(out.y+out.gam*Ax)/mugam, 1/mugam);
%         [out.z, ~] = prob.callg(out.y/out.gam + Ax, 1/out.gam);
    end
    out.iterations = dualout.iterations;
    if isfield(dualout, 'operations'), out.operations = dualout.operations; end
    out.residual = dualout.residual;
    out.ts = dualout.ts;
    out.prob = prob;
    out.preprocess = dualout.preprocess;
end
