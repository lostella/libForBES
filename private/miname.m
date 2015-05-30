% Copyright (C) 2015, Lorenzo Stella and Panagiotis Patrinos
%
% This file is part of ForBES.
% 
% ForBES is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% ForBES is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU Lesser General Public License for more details.
% 
% You should have received a copy of the GNU Lesser General Public License
% along with ForBES. If not, see <http://www.gnu.org/licenses/>.

function out = miname(prob, opt)
    t0 = tic();
    
    if nargin < 2, opt = []; end
    opt = ProcessOptions(opt);
    
    if nargin < 1, error('the problem structure must be provided as first argument'); end
    if ~isfield(prob, 'processed') || ~prob.processed, prob = ProcessSeparableProblem(prob, opt); end
    
    preprocess = toc(t0);
    
    if nargin > 1
        dualout = minfbe(dualprob, opt);
    else
        dualout = minfbe(dualprob);
    end
    
    out = GetPrimalOutput(prob, dualprob, dualout);
    out.preprocess = preprocess + dualout.preprocess;
    if nargin > 1
        out.opt = opt;
    end
    out.dual = dualout;
end
