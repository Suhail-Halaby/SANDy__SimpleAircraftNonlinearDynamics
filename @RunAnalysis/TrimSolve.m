function [ctrlTrim, vTrim, resNorm, f_e, exitflag, A, B, As, Bs, Aa, Ba] = TrimSolve(obj, allocSurf, ...
                            u_e, theta_e, ctrl_e, mixOpts, solveOpts, linOpts, satOpts)
arguments
        obj
        allocSurf (1,:) double {mustBeInteger, mustBePositive}
        u_e     (1,1) double
        theta_e (1,1) double
        ctrl_e  (:,1) double
        mixOpts.M   (:,:) double = []          % (#allocated) x (#virtual) mixing matrix
        solveOpts.tol   (1,1) double = 1e-10
        solveOpts.maxIt (1,1) double = 100
        linOpts.psi    (1,1) double = 0
        linOpts.height (1,1) double = 100
        linOpts.rho    (1,1) double = 1.225
        linOpts.g      (1,1) double = 9.81
        linOpts.gravInForce   (1,1) logical = true
        linOpts.momentAboutCG (1,1) logical = true
        satOpts.lbV (:,1) double = []          % lower bounds on VIRTUAL controls
        satOpts.ubV (:,1) double = []          % upper bounds on VIRTUAL controls
end
    linArgs = namedargs2cell(linOpts);
    n = numel(ctrl_e);
    a = allocSurf(:);
    alpha_e = 0;

    % Default M = identity => virtual controls ARE the physical allocated ones
    M = mixOpts.M;
    if isempty(M), M = eye(numel(a)); end
    nv = size(M, 2);                            % number of virtual controls

    lbV = satOpts.lbV;  if isempty(lbV), lbV = -inf(nv,1); end
    ubV = satOpts.ubV;  if isempty(ubV), ubV =  inf(nv,1); end

    % Recover virtual control from the initial physical guess: v0 = M \ c0
    c0phys = ctrl_e(a);
    v0 = M \ c0phys;
    v0 = min(max(v0, lbV), ubV);

    function [r, J] = trimResidual(v)
        ctrl    = ctrl_e;
        ctrl(a) = M * v;                        % virtual -> physical
        [~, ~, r] = obj.SystemJacobian(u_e, alpha_e, theta_e, ctrl, linArgs{:});
        if nargout > 1
            Bfull = obj.CtrlJacobian(u_e, alpha_e, theta_e, ctrl, linArgs{:});
            J     = Bfull(:, a) * M;            % chain rule: df/dv = (df/dc) M
        end
    end

    opts = optimoptions('lsqnonlin', 'Display','off', ...
        'SpecifyObjectiveGradient', true, ...
        'FunctionTolerance', solveOpts.tol, 'MaxIterations', solveOpts.maxIt);

    [vTrim, ~, f_e, exitflag] = lsqnonlin(@trimResidual, v0, lbV, ubV, opts);

    ctrlTrim    = ctrl_e;
    ctrlTrim(a) = M * vTrim;
    resNorm     = norm(f_e);

    atLimit = (vTrim <= lbV + eps) | (vTrim >= ubV - eps);
    if exitflag <= 0 || resNorm > sqrt(solveOpts.tol)
        if any(atLimit)
            warning("RunAnalysis:TrimSaturated", ...
                "Residual ||f|| = %.3g; %d of %d virtual controls saturated.", ...
                resNorm, nnz(atLimit), nv);
        else
            warning("RunAnalysis:TrimResidual", ...
                "Residual ||f|| = %.3g, no control saturated (exitflag %d).", resNorm, exitflag);
        end
    end

    [A, ~, f_e, As, Aa] = obj.SystemJacobian(u_e, alpha_e, theta_e, ctrlTrim, linArgs{:});
    [B, ~, ~,   Bs, Ba] = obj.CtrlJacobian (u_e, alpha_e, theta_e, ctrlTrim, linArgs{:});
    resNorm = norm(f_e);
end