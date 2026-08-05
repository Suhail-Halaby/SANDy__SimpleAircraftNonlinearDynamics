function [L,D,M,Cl,Cd,Cm] = ViternaCorrigan(~, alpha, delta, chord, surfChord, qStrip)
% simplePlate  Thin flat-plate section model, valid over the full +/-180 deg range.
%
% Attached branch : thin-airfoil theory      Cl = 2*pi*alpha
% Separated branch: Viterna-Corrigan / Hoerner normal-force plate
%                     CN = CDMAX*sin(alpha)
%                     Cl = CN*cos(alpha),  Cd = CN*sin(alpha)
% The two are joined by the Beard & McLain sigmoid blend.
%
% Refs:
%   Viterna & Corrigan, NASA CP-2230, 1982.
%   Hoerner, Fluid-Dynamic Lift, 1975 (CD ~ 2.0 for a 2-D normal plate).
%   Beard & McLain, Small Unmanned Aircraft, Princeton 2012, eqs 4.9-4.11.
%
% Inputs (angles in RADIANS):
%   alpha     local angle of attack                        [rad]
%   delta     control deflection, +ve = TE down            [rad]
%   chord     local chord c                                [m]
%   surfChord control-surface chord c_f (0 if none)        [m]
%   qStrip    local dynamic pressure 0.5*rho*V^2           [Pa]
%
% Outputs are PER UNIT SPAN (caller multiplies by strip width ds):
%   L,D,M     lift, drag, pitching moment about the 1/4-chord
%   Cl,Cd,Cm  section coefficients

CDMAX  = 2.0;       % normal-force plateau, 2-D plate broadside
CD0    = 0.02;      % skin friction at zero incidence
ASTALL = deg2rad(12);   % blend centre
MBLEND = 40;            % blend sharpness (larger = sharper stall)

% --- flap: thin-airfoil effectiveness, folded in as an incidence shift ---
E   = min(max(surfChord./chord, 0), 1);
th  = acos(2*E - 1);
Cl_delta = (surfChord > 0) .* (2*(pi - th + sin(th)));
Cm_delta = (surfChord > 0) .* (0.5*sin(th).*(cos(th) - 1));
tau = Cl_delta ./ (2*pi);

% --- wrap alpha to [-pi, pi] so the model is single-valued ---
a = atan2(sin(alpha), cos(alpha));

% --- separation blend: sigma = 0 attached, 1 fully separated ---
ep = exp(-MBLEND*(abs(a) - ASTALL));
sigma = 1 ./ (1 + ep);

% --- effective incidence (flap only bites while attached) ---
aEff = a + (1 - sigma).*tau.*delta;

% --- attached branch ---
Cl_att = 2*pi*aEff;
Cd_att = CD0;

% --- separated branch: CN = CDMAX*sin(a), resolved to wind axes ---
CN     = CDMAX * sin(a);
Cl_sep = CN .* cos(a);
Cd_sep = CN .* sin(a) + CD0;

% --- blended section coefficients ---
Cl = (1 - sigma).*Cl_att + sigma.*Cl_sep;
Cd = (1 - sigma).*Cd_att + sigma.*Cd_sep;

% --- moment: flap term (attached) + CP migration c/4 -> c/2 (separated) ---
CNtot = Cl.*cos(a) + Cd.*sin(a);
xcp   = 0.25 + 0.25*min(abs(a), pi - abs(a))/(pi/2);   % chords aft of LE
Cm    = (1 - sigma).*(Cm_delta.*delta) - sigma.*CNtot.*(xcp - 0.25);

% --- dimensional, per unit span ---
L = qStrip .* chord    .* Cl;
D = qStrip .* chord    .* Cd;
M = qStrip .* chord.^2 .* Cm;
end