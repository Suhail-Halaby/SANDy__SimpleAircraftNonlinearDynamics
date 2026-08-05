function [TotalF, TotalM, Mass, CG, Iv] = runVehicle(plane, ctrl, bodyRates, bodyVel, orQuat, height, rho, ModelName)
%RUNVEHICLE  One control-loop step: tilt -> inertia -> prop -> aero -> forcing.
%
%   [TotalF, TotalM, Mass, CG, Iv] = runVehicle(plane, ctrl, bodyRates, bodyVel, rho, ModelName)
%
%   plane - struct built once after vehicle config:
%       .wing, .wingL, .Vstab, .Hstab, .HstabL, .Fusl   StripModel objects
%       .propParams                                     base propellerBEMT params
%                                                         (D, P, B, rho, spin - NOT axis/hubloc/V_inf,
%                                                          those get refreshed every call below)
%       .wingWakeDist                                    distance from prop disk to wing, e.g. 0.25*D
%       .rpmFromThrottle                                 function handle: rpm = f(throttle)
%       .addedInertia                                    optional struct (Mass, MassLoc, Iv), or [] to skip
%
%   ctrl      
%   bodyRates = [p; q; r]
%   bodyVel   = [u; v; w]      <- body velocity. NOT the control vector.
%
%   Returns total body-frame force/moment, plus the Mass/CG/Iv computed
%   this step (tilting moves them, so the caller needs the current
%   values too, not just the ones from config time).

ailR   = ctrl(1);
ailL   = ctrl(2);
elv = ctrl(3);
throttle = ctrl(4);
throttleL = ctrl(5);
rudder = ctrl(6);
tiltAngle = ctrl(7);

if throttle > 1 
    warning('Throttle (R) Saturated at Max')
    throttle = 1;
end

if throttleL > 1 
    warning('Throttle (L) Saturated at Max')
    throttleL = 1;
end


%% 1-2. Tilt the surfaces that actually tilt

plane.wing.tiltSurf(tiltAngle);
plane.wingL.tiltSurf(tiltAngle);
% Vstab / Hstab / HstabL / Fusl don't tilt in this vehicle - add them
% here too if that ever changes.

%% 3. Re-evaluate inertia (tilting moves mass, so CG/Iv move with it)

[Mass, CG, Iv] = assembleInertia(plane.wing, plane.wingL, plane.Vstab, ...
    plane.Hstab, plane.HstabL, plane.Fusl);

if isfield(plane, 'addedInertia') && ~isempty(plane.addedInertia)
    ai  = plane.addedInertia;
    CG  = (Mass * CG + ai.Mass * ai.MassLoc) / (Mass + ai.Mass);
    Mass = Mass + ai.Mass;
    Iv   = Iv + ai.Iv;
end

%% 4. Apply prop thrust (mount geometry + inflow both change with tilt/speed)

pp        = plane.propParams;
pp.axis   = plane.wing.BoundProps.Axis;   % tilts with the wing
pp.hubloc = plane.wing.BoundProps.Hub;
pp.V_inf  = dot(bodyVel, pp.axis / norm(pp.axis));   % flow along prop axis

rpm  = plane.rpmFromThrottle(throttle);    % TODO: your "throttle forcing" item
prop = propellerBEMT(rpm, pp);
wake = prop.develop(plane.wingWakeDist);


% duplicate for other prop

ppL        = plane.propParams;
ppL.axis   = plane.wingL.BoundProps.Axis;   % tilts with the wing
ppL.hubloc = plane.wingL.BoundProps.Hub;
ppL.V_inf  = dot(bodyVel, ppL.axis / norm(ppL.axis));   % flow along prop axis
ppL.spin = -pp.spin;

rpmL  = plane.rpmFromThrottle(throttleL);    % TODO: your "throttle forcing" item
propL = propellerBEMT(rpmL, ppL);
wakeL = propL.develop(plane.wingWakeDist);



%% 5. Aerodynamics -> total forcing

[F,  M ] = plane.wing.aeromodel (CG, ailR, bodyVel, bodyRates, rho, ModelName, ...
    wake.va', wake.vt', wake.r' ./ wake.R);

[FL, ML] = plane.wingL.aeromodel(CG, ailL, bodyVel, bodyRates, rho, ModelName, ...
    wakeL.va', wakeL.vt', wakeL.r' ./ wakeL.R);

[Fv, Mv] = plane.Vstab.aeromodel(CG, rudder,  bodyVel, bodyRates, rho, ModelName);

[Fh, Mh] = plane.Hstab.aeromodel(CG, elv,  bodyVel, bodyRates, rho, ModelName);

[Fhl, Mhl] = plane.HstabL.aeromodel(CG, elv,  bodyVel, bodyRates, rho, ModelName);


% gravitation contribution:
orQuat = quaternion(orQuat(1), orQuat(2), orQuat(3), orQuat(4));

GravF = rotmat(orQuat,'frame')*[0;0;Mass*9.81];


TotalF = GravF + F + FL + Fv + Fh + Fhl + prop.T_axis + propL.T_axis;
TotalM =         M + ML + Mv + Mh + Mhl + prop.M_axis + propL.M_axis;

[TotalF, TotalM] = applyFloorContact(TotalF, TotalM, orQuat, bodyVel, height);

end