function obj = tiltSurf(obj, tilt)
% tiltSurf  Set the Tilt angle [deg], re-discretise, and re-apply everything
%           that was derived from the strip grid.
%
%   wing.tiltSurf(60)
%
% Tilt rotates the panel about body Y through RootQC+TiltOffset. The strip
% grid is rebuilt at the same resolution/spacing, then:
%   - bound props are replayed (Mask/Cover/Swirl are sized to the old grid,
%     and surface-mounted props rotate with the wing)
%   - the mass model is replayed (strip areas/mid-chords move, so PartCG and
%     Inertias change - this is the inertial contribution of the tilt)
%
% Props with AxisFixed=true, and point masses with LocFixed=true, keep their
% body-axis definition and do not follow the tilt.
%
% Both refresh calls are no-ops if the corresponding add* was never called.
arguments
    obj
    tilt (1,1) double
end

obj.Tilt = tilt;

if isfield(obj.Strips,'n')
    obj.discretise('nStrips', obj.Strips.n, 'Spacing', obj.Strips.spacing);
else
    obj.discretise();
end

obj.refreshBoundProps();   % rebuild props on the new grid / tilted frame
obj.refreshMass();         % rebuild mass, CG, inertia
end