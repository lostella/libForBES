function [fz, gradfz] = LogReg(z)
    m = length(z);
    pz = 1./(1+exp(-z));
    fz = -sum(log(pz))/m;
    if nargout >= 2
        gradfz = (pz-1)/m;
    end
end
