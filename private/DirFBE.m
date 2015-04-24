function [cachet, cnt] = DirFBE(prob, gam, tau, cache, mode, cachet)
% Computes the 'directional' FBE, i.e., FBE(x+tau*d) and its derivative with
% respect to tau, if requested. Here x = cache.x, d = cache.dir.
% If cachet (6th argument) is provided, then skips precomputing the data
% that has already been stored in cachet.
%
% If mode == 1, then compute only FBE(x+tau*d) and put it into cachet.
% If mode == 2, compute only dFBE(x+tau*d), the directional derivative.
% If mode == 3, compute both FBE and dFBE at x+tau*d.
    
    %      Q, A, C, f2, g
    cnt = [0, 0, 0, 0, 0];

    if nargin < 6
        fxt = 0;
        gradfxt = 0;
        if prob.istheref1
            cachet.res1x = cache.res1x + tau*cache.Adir;
            cachet.Qres1x = cache.Qres1x + tau*cache.QAdir;
            cachet.gradf1x = cache.gradf1x + tau*cache.ATQAdir;
            cachet.f1x = cache.f1x + tau*cache.f1linear + (0.5*tau^2)*cache.f1quad;
            fxt = fxt + cachet.f1x;
            gradfxt = gradfxt + cachet.gradf1x;
        end
        if prob.istheref2
            cachet.res2x = cache.res2x + tau*cache.Cdir;
            if prob.useHessian
                [f2xt, gradf2res2xt, cachet.Hessf2res2x] = prob.f2(cachet.res2x);
            else
                [f2xt, gradf2res2xt] = prob.f2(cachet.res2x);
            end
            cnt(4) = cnt(4)+1;
            if prob.isthereC
                if prob.isCfun, gradf2xt = prob.CT(gradf2res2xt);
                else gradf2xt = prob.C'*gradf2res2xt; end
                cnt(3) = cnt(3)+1;
            else
                gradf2xt = gradf2res2xt;
            end
            fxt = fxt + f2xt;
            gradfxt = gradfxt + gradf2xt;
        end
        if prob.istherelin
            cachet.flinx = cache.flinx + tau*cache.lindir;
            fxt = fxt + cachet.flinx;
            gradfxt = gradfxt + prob.l;
        end
        % compute proximal gradient step
        cachet.x = cache.x + tau*cache.dir;
        cachet.fx = fxt;
        cachet.gradfx = gradfxt;
        yt = cachet.x - gam*gradfxt;
        [cachet.z, cachet.gz] = prob.g(yt, gam);
        cnt(5) = cnt(5)+1;
        cachet.diff = cachet.z-cachet.x;
    end
    
    if mode == 1 || mode == 3
        sqnormdifft = cachet.diff'*cachet.diff;
        cachet.normdiff = sqrt(sqnormdifft);
        cachet.FBE = cachet.fx + cachet.gz + ...
            cachet.gradfx'*cachet.diff + ...
            (0.5/gam)*sqnormdifft;
    end
    
    if mode >= 2
        Hdir = 0;
        if prob.istheref1
            Hdir = Hdir + cache.ATQAdir;
        end
        if prob.istheref2
            if prob.useHessian
                HCdir = cachet.Hessf2res2x(cache.Cdir);
            else
                res2xtepsdir = cachet.res2x + 1e-100i*cache.Cdir;
                [~, gradf2res2xtepsdir] = prob.f2(res2xtepsdir);
                cnt(4) = cnt(4)+1;
                HCdir = imag(gradf2res2xtepsdir)/1e-100;
            end
            if prob.isthereC
                if prob.isCfun, Hdir = Hdir + prob.CTfun(HCdir);
                else Hdir = Hdir + (prob.C'*HCdir); end
                cnt(3) = cnt(3)+1;
            else
                Hdir = Hdir + HCdir;
            end
        end
        cachet.dFBE = cachet.diff'*Hdir-(cachet.diff'*cache.dir)/gam;
    end
end
