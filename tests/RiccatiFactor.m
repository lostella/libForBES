function [LRs, Ks, Ms, Ls] = RiccatiFactor(Q, R, Qf, A, B, N)
    n = size(Q,1);
    m = size(R,1);
    Ps = zeros(n, n, N+1);
    Ps(:,:,N+1) = Qf;
    LRs = zeros(m, m, N);
    Ss = zeros(m, n, N);
    Ks = zeros(m, n, N);
    Ms = zeros(m, n, N);
    Ls = zeros(n, n, N);
    for k = N:-1:1
        Rbar = R+B'*(Ps(:,:,k+1)*B);
        LR = chol(Rbar, 'lower');
        LRs(:,:,k) = LR;
        Ss(:,:,k) = B'*(Ps(:,:,k+1)*A);
        Ks(:,:,k) = -(LR'\(LR\Ss(:,:,k)));
        Ps(:,:,k) = Q + A'*(Ps(:,:,k+1)*A) + Ss(:,:,k)'*Ks(:,:,k);
    end
    for k = 1:N
        LR = LRs(:,:,k);
        Ms(:,:,k) = -(LR'\(LR\B'));
        Ls(:,:,k) = (A+B*Ks(:,:,k))';
    end
end
