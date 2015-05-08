function obj = separableSum(idx, sizes, objs)
    obj.makeprox = @() make_separableSum_prox(idx, sizes, objs);
end

function op = make_separableSum_prox(idx, sizes, objs)
    proxes = {};
    for i=1:length(objs)
        proxes{end+1} = objs{i}.makeprox();
    end
    op = @(x, gam) call_separableSum_prox(x, gam, idx, sizes, proxes);
end

function [prox, val] = call_separableSum_prox(x, gam, idx, sizes, proxes)
    n = sum(sizes);
    prox = zeros(n, 1);
    val = 0;
    baseidx = 0;
    for i=1:length(idx)
        endcurr = baseidx+sizes(i);
        xcurr = x(baseidx+1:endcurr);
        [prox(baseidx+1:endcurr), val1] = proxes{idx(i)}(xcurr, gam);
        val = val+val1;
        baseidx = endcurr;
    end
end
