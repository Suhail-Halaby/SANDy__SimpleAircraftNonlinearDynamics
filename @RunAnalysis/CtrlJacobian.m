function [B, x_e, f_e, Bs, Ba] = CtrlJacobian(obj, u_e, alpha_e, theta_e, ctrl_e, opts)
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
nc = size(ctrl_e,1);
n = 9;
B = zeros(n,nc);


for j = 1:nc
    e  = zeros(nc,1);
    e(j) = 1;
    hj = 0.005;
    fpp = obj.stateDeriv(plane, x_e, ctrl_e  + 2*hj*e, cfg);
    fp  = obj.stateDeriv(plane, x_e, ctrl_e  +   hj*e, cfg);
    fm  = obj.stateDeriv(plane, x_e, ctrl_e  -   hj*e, cfg);
    fmm = obj.stateDeriv(plane, x_e, ctrl_e  - 2*hj*e, cfg);
    B(:,j) = (-fpp + 8*fp - 8*fm + fmm) / (12*hj);
end


% partion into symmetric and asymmetric:
symStates = [1;3;5;8]; % u w q theta
asymStates = [2;4;6;7;9]; % v p r phi psi

Bs = B(symStates,:);
Ba = B(asymStates,:);
end