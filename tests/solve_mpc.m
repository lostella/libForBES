function out = solve_mpc(Q, R, Qf, A, B, F, G, Ff, c, cf, weights, N, x0, opt, solver)
%SOLVE_MPC Sets up the QP associated with the specified MPC problem and
%   solves it using MINAME/MINFBE. The arguments refer to the following MPC
%   problem formulation:
%
%       minimize sum{ 0.5 x_i'Q x_i + 0.5 u_i'R u_i, i=0,...,N-1 }
%                    + 0.5 x_N'Q_f x_N + g(z, z_N)
%
%       subject to  x_{i+1} = A x_i + B u_i,    i = 0,...,N-1
%                   F x_i + G u_i + z = c,      i = 0,...,N-1
%                   F_f x_N + z_N = c_f
%                   x_0 = x0 (given)
%
%   Here variables z, z_N have the role of slack variables, and they are
%   penalized in the g term of the cost through some weights:
%
%       g(z,z_N) = -w.*max{0,z}-w_N.*max{0,z_N}
%
%   and w, w_N are provided in the weights = [w; w_N] vector. If weights is
%   a scalar s then it is assumed that w = w_N = [s; ...; s], i.e., all weights
%   are the same. If some weight equals +inf then the corresponding constraint
%   is a (hard) nonnegativity constraint.
%
    t0 = tic;
    
    if nargin < 14
        opt = [];
    end
    
    % this is just in case one wants to plug some other solver
    if nargin < 15
        solver = @(p,o) miname(p,o);
    end
    
    [n, m] = size(B);
    
    % make big (sparse) constraint matrix blocks
    FG = sparse([F, G]);
    constrBlocks = cell(N+1,1);
    for i=1:N
        constrBlocks{i} = FG;
    end
    constrBlocks{N+1} = sparse(Ff);
    
    % pack up problem
    prob.A1 = blkdiag(constrBlocks{:});
    [LRs, Ks, Ms, Ls] = RiccatiFactor(Q, R, Qf, A, B, N);
    prob.x1step = @(w) RiccatiSolve(w, x0, A, B, LRs, Ks, Ms, Ls, int32(n), int32(m), int32(N));
    prob.B = 1;
    prob.c = [repmat(c, N, 1); cf];
    prob.zstep = @(y, gam) SoftNonnegativeIndicator(-y/gam, weights, 1/gam);
    
    preprocess = toc(t0);
    
    if isempty(opt)
        primalout = solver(prob);
    else
        primalout = solver(prob, opt);
    end
    
    out.x = primalout.x1;
    out.iterations = primalout.iterations;
    out.time = primalout.ts(end);
    out.preprocess = preprocess + primalout.preprocess;
    out.primal = primalout;
end
