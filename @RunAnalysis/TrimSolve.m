function [ctrlTrim, A, B, resNorm, f_e, As, Bs, Aa, Ba] = TrimSolve(obj, allocSurf, ...
                            u_e, alpha_e, theta_e, ctrl_e, solveOpts, linOpts)
    arguments
        obj
        allocSurf (1,:) double {mustBeInteger, mustBePositive}
        u_e     (1,1) double 
        alpha_e (1,1) double
        theta_e (1,1) double
        ctrl_e  (:,1) double
        solveOpts.tol   (1,1) double = 1e-4
        solveOpts.maxIt (1,1) double = 50
        linOpts.psi    (1,1) double = 0
        linOpts.height (1,1) double = 100
        linOpts.rho    (1,1) double = 1.225
        linOpts.g      (1,1) double = 9.81
        linOpts.gravInForce   (1,1) logical = true
        linOpts.momentAboutCG (1,1) logical = true
    end

    linArgs = namedargs2cell(linOpts);   % expand struct -> name-value for the Jacobians

    ctrlTrim  = ctrl_e;                   % current control guess
    resNorm   = inf;
    converged = false;

    for it = 1:solveOpts.maxIt
        % Residual and control effectiveness at the CURRENT controls
        [~, ~, f_e] = obj.SystemJacobian(u_e, alpha_e, theta_e, ctrlTrim, linArgs{:});
        B           = obj.CtrlJacobian (u_e, alpha_e, theta_e, ctrlTrim, linArgs{:});

        resNorm = norm(f_e);
        if resNorm <= solveOpts.tol
            converged = true;
            break;
        end

        % Newton step on the allocated surfaces only:  f_e + B_trim*du = 0
        B_trim        = B(:, allocSurf);
        du            = zeros(numel(ctrlTrim), 1);
        du(allocSurf) = -pinv(B_trim) * f_e;

        ctrlTrim = ctrlTrim + du;
    end

    if ~converged
        warning("RunAnalysis:TrimNotConverged", ...
            "Trim did not converge: ||f|| = %.3g after %d iters (tol %.3g).", ...
            resNorm, solveOpts.maxIt, solveOpts.tol);
    end

    % Final linear model AT the trim point (also refreshes resNorm honestly)
    [A, ~, f_e, As, Aa] = obj.SystemJacobian(u_e, alpha_e, theta_e, ctrlTrim, linArgs{:});
    [B, ~, ~,   Bs, Ba] = obj.CtrlJacobian (u_e, alpha_e, theta_e, ctrlTrim, linArgs{:});
    resNorm = norm(f_e);
end