function out = forbes(prob, opt)
    prob = IdentifyProblem(prob);
    switch prob.identified
        case 1
            if nargin > 1
                out = minfbe(prob, opt);
            else
                out = minfbe(prob);
            end
        case 2
            if nargin > 1
                out = miname(prob, opt);
            else
                out = miname(prob);
            end
    end
end

function prob = IdentifyProblem(prob)
    % Check whether we are given equality constraints or not
    if any(isfield(prob, {'A1', 'A1t', 'A2', 'A2t', 'B', 'b'}))
        flagEQ = true;
    else
        flagEQ = false;
    end
    % Check whether we are given affine mappings composed with f1 and f2
    if any(isfield(prob, {'C1', 'C1t', 'd1', 'C2', 'C2t', 'd2'}))
        flagAFF = true;
    else
        flagAFF = false;
    end
    % Check for uncertain situations
    if (flagEQ && flagAFF)
        error('you cannot provide both equality constraints and affine mappings composed with f1, f2');
    end
    if flagEQ
        % If equality constraints are provided, solve the dual
        prob.identified = 2;
    else
        % Otherwise solve the primal
        prob.identified = 1;
    end
end
