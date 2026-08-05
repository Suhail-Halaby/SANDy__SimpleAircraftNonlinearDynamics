function [prop] = propellerBEMT(rpm,params)
% PROPELLERBEMT  BEMT prop model (assumed geometric twist), valid to the
% static limit, with McCormick-style streamwise slipstream development so the
% induced-flow field can be sampled at any axial station behind the disk.
%
% INPUTS
%   rpm    : propeller speed [rev/min]
%   params : struct
%       .D, .P, .B, .rho, .nElem, .rHub, .Cla, .Cd0, .chord   (as before)
%       .axis    thrust/forward axis (unit or not), default [1;0;0]
%       .hubloc  hub position rel. to CG [m], default [0;0;0]
%       .spin    rotation sense about +axis: +1 (RH) or -1 (LH), default +1
%       .x_eval  OPTIONAL streamwise station(s) [m] to evaluate the wake at
%
% OUTPUT struct 'prop'
%   ... disk-plane fields (r, va, vt, a, aprime, T, Q, CT, CP, ...) ...
%   .T_axis  thrust force vector along axis            [3x1, N]
%   .Q_axis  reaction torque on airframe about axis    [3x1, N*m]
%   .M_hub   moment from thrust offset about CG         [3x1, N*m]
%   .M_axis  total moment about CG (M_hub + Q_axis)     [3x1, N*m]
%   .develop function handle: S = prop.develop(s)
%   .slipstream populated if params.x_eval given

if ~isfield(params,'rho'),    params.rho    = 1.225;    end
if ~isfield(params,'nElem'),  params.nElem  = 40;       end
if ~isfield(params,'rHub'),   params.rHub   = 0.15;     end
if ~isfield(params,'Cla'),    params.Cla    = 2*pi;     end
if ~isfield(params,'Cd0'),    params.Cd0    = 0.02;     end
if ~isfield(params,'axis'),   params.axis   = [1;0;0];  end
if ~isfield(params,'hubloc'), params.hubloc = [0;0;0];  end
if ~isfield(params,'spin'),   params.spin   = 1;        end

axis   = params.axis / norm(params.axis);
hubloc = params.hubloc;                 % <-- was missing; bodyVel path needed it
spin   = params.spin;

if ~isfield(params,'V_inf') && ( ~isfield(params,'bodyVel') &&  ~isfield(params,'bodyRates') )
    error("No Velocity Supplied By User")
elseif ~isfield(params,'V_inf') && ( isfield(params,'bodyVel') &&  isfield(params,'bodyRates') )
    warning("Using Axis-Resolved Vehicle Motion")
    vlocal = params.bodyVel + cross(params.bodyRates,hubloc);
    V_inf = dot(axis,vlocal);
elseif isfield(params,'V_inf') && ( ~isfield(params,'bodyVel') &&  ~isfield(params,'bodyRates') )
    V_inf = params.V_inf;
else
    warning("Both velocities supplied, defaulting to V_inf")
    V_inf = params.V_inf;
end

R     = params.D/2;
B     = params.B;
rho   = params.rho;
n     = rpm/60;
Omega = 2*pi*n;

if Omega <= 0
    warning('rpm must be > 0 , setting rpm to 1 rpm');
    Omega = 2*pi / 60;
end

rR = linspace(params.rHub, 0.99, params.nElem).';
r  = rR * R;

if isfield(params,'chord')
    if isa(params.chord,'function_handle'), c = params.chord(rR);
    else,                                   c = params.chord*ones(size(rR));
    end
else
    c = (0.24 - 0.10*rR) * R;
end

beta = atan2(params.P, 2*pi*r);

a      = zeros(size(r));
aprime = zeros(size(r));
va     = zeros(size(r));   % axial induced velocity at disk
dT     = zeros(size(r));
dQ     = zeros(size(r));

maxIt = 300; tol = 1e-6; relax = 0.25;
STATIC_TOL = 1e-3;

for i = 1:numel(r)
    if abs(V_inf) < STATIC_TOL
        % ---------- STATIC / HOVER BRANCH ----------
        vai = 0.5*Omega*r(i)*0.1 + 1;
        api = 0.01;
        for it = 1:maxIt
            Ua = vai;
            Ut = Omega*r(i)*(1 - api);
            phi = atan2(Ua, Ut);
            W2  = Ua^2 + Ut^2;

            alpha = beta(i) - phi;
            Cl = params.Cla*alpha;
            Cd = params.Cd0 + 0.02*alpha^2;

            f = (B/2)*(1 - rR(i))/max(sin(phi),1e-3);
            F = max((2/pi)*acos(min(1,exp(-f))), 1e-3);

            cphi = cos(phi); sphi = sin(phi);
            dTi = 0.5*rho*W2*B*c(i)*(Cl*cphi - Cd*sphi);
            dQi = 0.5*rho*W2*B*c(i)*(Cl*sphi + Cd*cphi)*r(i);

            va_new = sqrt(max(dTi,0) / (4*pi*r(i)*rho*F + eps));
            ap_new = dQi / (4*pi*r(i)^3*rho*Omega*max(vai,eps)*F + eps);

            va_new = max(min(va_new, 2*Omega*r(i)), 0);
            ap_new = max(min(ap_new, 0.9), 0);

            if abs(va_new-vai)<tol && abs(ap_new-api)<tol
                vai = va_new; api = ap_new; break;
            end
            vai = (1-relax)*vai + relax*va_new;
            api = (1-relax)*api + relax*ap_new;
        end
        va(i) = vai; aprime(i) = api;
        a(i)  = NaN;
        dT(i) = dTi; dQ(i) = dQi;
    else
        % ---------- ADVANCING BRANCH ----------
        ai = 0.1; api = 0.01;
        for it = 1:maxIt
            Ua = V_inf*(1 + ai);
            Ut = Omega*r(i)*(1 - api);
            phi = atan2(Ua, Ut);
            W2  = Ua^2 + Ut^2;

            alpha = beta(i) - phi;
            Cl = params.Cla*alpha;
            Cd = params.Cd0 + 0.02*alpha^2;

            f = (B/2)*(1 - rR(i))/max(sin(phi),1e-3);
            F = max((2/pi)*acos(min(1,exp(-f))), 1e-3);

            cphi = cos(phi); sphi = sin(phi);
            dTi = 0.5*rho*W2*B*c(i)*(Cl*cphi - Cd*sphi);
            dQi = 0.5*rho*W2*B*c(i)*(Cl*sphi + Cd*cphi)*r(i);

            a_new  = dTi / (4*pi*r(i)*rho*V_inf^2*(1+ai)*F + eps);
            ap_new = dQi / (4*pi*r(i)^3*rho*Omega*V_inf*(1+ai)*F + eps);

            a_new  = max(min(a_new, 1.5), -0.5);
            ap_new = max(min(ap_new, 0.9), -0.5);

            if abs(a_new-ai)<tol && abs(ap_new-api)<tol
                ai = a_new; api = ap_new; break;
            end
            ai  = (1-relax)*ai  + relax*a_new;
            api = (1-relax)*api + relax*ap_new;
        end
        a(i) = ai; aprime(i) = api;
        va(i) = ai*V_inf;
        dT(i) = dTi; dQ(i) = dQi;
    end
end

% ---- tangential (swirl) induced velocity at disk ----
vt = aprime .* Omega .* r;

% ---- integrate totals ----
T = trapz(r, dT);
Q = trapz(r, dQ);
P_shaft = Q*Omega;

% ---- resolve loads into the axis frame ----
T_axis = T * axis;                 % thrust force along +axis
Q_axis = -spin * Q * axis;         % airframe reaction torque (opposes prop spin)
M_hub  = cross(hubloc, T_axis);    % moment from thrust acting off the CG
M_axis = M_hub + Q_axis;           % total moment about CG

CT = T / (rho*n^2*params.D^4);
CP = P_shaft / (rho*n^3*params.D^5);

% ---- thrust-weighted disk averages ----
w = max(dT,0);
if trapz(r,w) > 0
    va_mean = trapz(r, va.*w)/trapz(r,w);
    vt_mean = trapz(r, vt.*w)/trapz(r,w);
else
    va_mean = mean(va); vt_mean = mean(vt);
end

% ---- fully-developed far wake (s -> inf limit, kept for compatibility) ----
vs        = V_inf + 2*va;
Rcontract = R * sqrt( max(V_inf + va_mean,eps) / max(V_inf + 2*va_mean,eps) );

% ================= McCormick streamwise slipstream development =================
% Vortex-tube development factor grows from the disk value (s=0) to twice the
% disk value (s->inf):   f(s) = 1 + s/sqrt(s^2 + R^2)
    function S = slipstreamAt(s)
        s  = s(:).';
        fs = 1 + s ./ sqrt(s.^2 + R^2);

        S.s      = s;
        S.f      = fs;
        S.va     = va * fs;
        S.Vaxial = V_inf + va * fs;
        S.vt     = vt * fs;
        S.r      = r .* sqrt( (V_inf + va) ./ (V_inf + va.*fs + eps) );
        S.R      = R * sqrt( (V_inf + va_mean) ./ (V_inf + va_mean.*fs + eps) );
    end

% ---- pack output ----
prop = struct('r',r,'rR',rR,'beta',beta,'chord',c, ...
              'a',a,'aprime',aprime,'va',va,'vt',vt, ...
              'vs',vs,'Rcontract',Rcontract, ...
              'dT',dT,'dQ',dQ,'T',T,'Q',Q,'P_shaft',P_shaft, ...
              'CT',CT,'CP',CP,'va_mean',va_mean,'vt_mean',vt_mean, ...
              'Omega',Omega,'n',n,'V_inf',V_inf, ...
              'T_axis',T_axis,'Q_axis',Q_axis,'M_hub',M_hub,'M_axis',M_axis,'spin',spin);

prop.develop = @slipstreamAt;
if isfield(params,'x_eval') && ~isempty(params.x_eval)
    prop.slipstream = slipstreamAt(params.x_eval);
end
end