function [TotalF, TotalM] = applyFloorContact(TotalF, TotalM, orQuat, bodyVel, height)
%APPLYFLOORCONTACT  Trivial rigid-floor constraint.
%
%   If height <= 0, cancel any net body-frame force that would push the
%   vehicle further downward (i.e. accelerate it below the floor), and
%   also catch the case where it's still carrying downward velocity from
%   the previous step. This is NOT a real contact/landing-gear model -
%   no friction, no bounce, no moment reaction - just "can't sink below
%   the floor."
%
%   Convention: assumes world frame is NED-style (+Z DOWN), matching the
%   GravF = rotmat(orQuat,'frame')*[0;0;Mass*9.81] line in runVehicle.
%   "height" should be passed in as altitude, positive UP. If your world
%   frame is +Z UP instead, flip the two sign checks marked below.
groundTol = 1e-4;   % meters
if height > groundTol
    return;
end

% rotmat(orQuat,'frame') rotates world -> body (used for GravF above).
% body -> world is just the transpose (rotation matrices are orthogonal).
R_wb = rotmat(orQuat, 'frame');   % world -> body
R_bw = R_wb.';                    % body  -> world

F_world   = R_bw * TotalF;
vDownWorld = R_bw(3,:) * bodyVel;   % world-frame downward velocity

% +Z is DOWN in world frame here, so "into the floor" = positive Z
if F_world(3) > 0 || vDownWorld > 0
    F_world(3) = 0;                % kill the downward force component
    TotalF = R_wb * F_world;       % rotate back to body frame
end

% No moment correction here (trivial model) - TotalM passed through.


end