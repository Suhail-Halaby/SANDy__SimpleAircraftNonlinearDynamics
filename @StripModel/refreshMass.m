function [] = refreshMass(obj)
% refreshMass  Rebuild the mass/CG/inertia model against the current strip
%              grid and the current Tilt.
%
%   obj.refreshMass()
%
%   discretise() rebuilds obj.Strips, which invalidates PartMass, PartCG and
%   Inertias: strip areas and mid-chord positions change with the grid, and
%   every position was computed in the OLD tilted frame. This replays the
%   stored MassSpec through addMass so the inertia follows the surface.
%
%   Call after ANY re-discretise - tiltSurf() already does.
%   No-op if addMass has never been called.

if ~isprop(obj,'MassSpec') || isempty(obj.MassSpec)
    return;
end

if isempty(obj.Strips) || ~isfield(obj.Strips,'n')
    error('refreshMass:NoStrips', ...
        'obj.Strips is empty - call discretise() before refreshing mass.');
end

% addMass overwrites rather than appends, so no clear step is needed here
% (unlike refreshBoundProps, where replaying would duplicate props).
obj.addMass(obj.MassSpec{:});
end