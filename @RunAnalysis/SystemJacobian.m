function [A, x_e, f_e, As, Aa] = SystemJacobian(obj, u_e, alpha_e, theta_e, ctrl_e, opts)
arguments
    obj
    u_e     (1,1) double 
    alpha_e (1,1) double                       % equilibrium angle of incidence
    theta_e (1,1) double                       % equilibrium pitch [rad]
    ctrl_e  (:,1) double                       % equilibrium controls
    opts.psi    (1,1) double = 0
    opts.height (1,1) double = 100
    opts.rho    (1,1) double = 1.225
    opts.g      (1,1) double = 9.81
    opts.gravInForce    (1,1) logical = true    % true = runVehicle's F already includes weight
    opts.momentAboutCG  (1,1) logical = true    % true = runVehicle's M is taken about the CG
    opts.Steps  (1,9) double = [0.1 0.1 0.1  0.05 0.05 0.05  0.005 0.005 0.005]
end

plane = obj.loadVehicle();
assert(isfield(plane,'ModelName'), "plane.ModelName missing after buildVehicle");

cfg = struct('height',opts.height,'rho',opts.rho,'g',opts.g, ...
    'gravInForce',opts.gravInForce,'momentAboutCG',opts.momentAboutCG);

% --- Equilibrium state: wings-level, no sideslip, alpha = theta - gamma ---
w_e     = u_e * tan(alpha_e);
x_e = [ u_e;   % u
    0;                  % v
    w_e;   % w
    0; 0; 0;            % p q r
    0;                  % phi
    theta_e;            % theta
    opts.psi ];         % psi

% Trim residual — should be ~0 if x_e is a genuine equilibrium
f_e = obj.stateDeriv(plane, x_e, ctrl_e, cfg);

% --- 5-point central difference:  f'(x) = (-f(x+2h)+8f(x+h)-8f(x-h)+f(x-2h))/(12h) ---
n = 9;
A = zeros(n, n);
h = opts.Steps(:);
for j = 1:n
    e  = zeros(n,1); e(j) = 1;
    hj = h(j);
    fpp = obj.stateDeriv(plane, x_e + 2*hj*e, ctrl_e, cfg);
    fp  = obj.stateDeriv(plane, x_e +   hj*e, ctrl_e, cfg);
    fm  = obj.stateDeriv(plane, x_e -   hj*e, ctrl_e, cfg);
    fmm = obj.stateDeriv(plane, x_e - 2*hj*e, ctrl_e, cfg);
    A(:,j) = (-fpp + 8*fp - 8*fm + fmm) / (12*hj);
end


% partion into symmetric and asymmetric:
symStates = [1;3;5;8]; % u w q theta
asymStates = [2;4;6;7;9]; % v p r phi psi

As = A(symStates,symStates);
Aa = A(asymStates,asymStates);
end