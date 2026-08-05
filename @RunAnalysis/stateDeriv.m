function xdot = stateDeriv(obj, plane, x, ctrl, cfg)
% State x = [u v w p q r phi theta psi]', body-frame vels/rates + ZYX Euler
u=x(1); v=x(2); w=x(3);
p=x(4); q=x(5); r=x(6);
phi=x(7); theta=x(8); psi=x(9); 

bodyVel   = [u; v; w];
bodyRates = [p; q; r];
orQuat    = obj.euler2quat(phi, theta, psi);

[F, M, Mass, CG, Iv] = obj.RunFcn(plane, ctrl, bodyRates, bodyVel, ...
    orQuat, cfg.height, cfg.rho, plane.ModelName);
F = F(:); M = M(:); CG = CG(:);

% Add weight only if runVehicle's force does NOT already include it
if ~cfg.gravInForce
    g_body = cfg.g * [ -sin(theta); cos(theta)*sin(phi); cos(theta)*cos(phi) ];
    F = F + Mass * g_body;
end

% Moments must act about the CG
if ~cfg.momentAboutCG
    M = M - cross(CG, F);
end

% Newton–Euler in body frame
uvw_dot = F ./ Mass - cross(bodyRates, bodyVel);
pqr_dot = Iv \ (M - cross(bodyRates, Iv * bodyRates));

% ZYX Euler kinematics
sphi = sin(phi); cphi = cos(phi);
cth  = cos(theta); tth = tan(theta);
phi_dot   = p + tth * (q*sphi + r*cphi);
theta_dot =      q*cphi - r*sphi;
psi_dot   =     (q*sphi + r*cphi) / cth;

xdot = [uvw_dot; pqr_dot; phi_dot; theta_dot; psi_dot];
end