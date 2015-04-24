function [cachet, t, cnt, exitflag] = ArmijoLS(prob, gam, cache, df0, lsopt, fr)
    %      Q, A, C, f2, proxg
    cnt = [0, 0, 0, 0, 0];
    t = lsopt.tau0;
    arm_hi = lsopt.delta*df0;
    exitflag = -1;
    if nargin >= 6
        f0 = fr;
    else
        f0 = cache.FBE;
    end
    for i = 1:lsopt.nLS
        [cachet, cnt1] = DirFBE(prob, gam, t, cache, 1);
        cnt = cnt+cnt1;
        ft = cachet.FBE;
        if ft <= f0 + t*arm_hi
            exitflag = 0;
            break;
        end
        if i == 1%quadratic interpolation
            tn = ArmijoQuadInterp(f0,df0,t,ft);
        else%cubic interpolation
            tn = ArmijoCubInterp(f0,df0,told,ftold,t,ft);
        end
        if tn <= 0
            tn = 0.5*t;
        end
        told = t;ftold = ft;
        t = tn;
        if t <= lsopt.progTol
            exitflag = -2;
            break
        end
    end
end

function t = ArmijoQuadInterp(f0,df0,t,ft)
    % Minimizer of interpolant belongs to [0,t1]
    tdf0 = t*df0;
    q = ft-f0-tdf0;
    if q > 0%quadratic is strongly convex
        t = -(tdf0*t)/(2*q);
    else
        t = -1;
    end
end

function t = ArmijoCubInterp(f,df,t0,f0,t1,f1)
    % Minimizer of interpolant belongs to [0,t1]
    t02 = t0^2;
    t12 = t1^2;
    ab = 1/(t02*t12*(t1-t0))*[t02 -t12;-t0^3 t1^3]*[f1-f-df*t1;f0-f-df*t0];
    a = ab(1);
    b = ab(2);
    t = (-b+sqrt(b^2-3*a*df))/(3*a);
end
