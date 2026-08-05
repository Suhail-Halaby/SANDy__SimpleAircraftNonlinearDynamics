function [L,D,M,Cl,Cd,Cm] = simplePlate(~, alpha, delta, chord, surfChord, qStrip)
% simplePlate  Thin-flat-plate section lifting model (private).
%
% Inputs (angles in RADIANS):
%   alpha     local angle of attack            [rad]
%   delta     control deflection, +ve = TE down (=> +lift)   [rad]
%   chord     local chord  c                   [m]
%   surfChord control-surface chord  c_f       [m]   (0 if no control here)
%   qStrip    local dynamic pressure 0.5*rho*V^2 [Pa]
%
% Outputs are PER UNIT SPAN (caller multiplies by strip width ds):
%   L,D,M     lift, drag, pitching moment about the 1/4-chord
%   Cl,Cd,Cm  section coefficients
    E  = min(max(surfChord./chord, 0), 1);
    th = acos(2*E - 1);

    Cl_delta = (surfChord > 0) .* (2*(pi - th + sin(th)));
    Cm_delta = (surfChord > 0) .* (0.5*sin(th).*(cos(th) - 1));

    % --- section coefficients (linear regime) ---
    Cl = 2*pi*alpha + Cl_delta.*delta;
    Cd = 0.0;                                      % inviscid plate, full LE suction
    Cm = Cm_delta.*delta;                          % flat plate => no alpha term at c/4

    % --- dimensional, per unit span ---
    L = qStrip .* chord    .* Cl;
    D = qStrip .* chord    .* Cd;
    M = qStrip .* chord.^2  .* Cm;
end
