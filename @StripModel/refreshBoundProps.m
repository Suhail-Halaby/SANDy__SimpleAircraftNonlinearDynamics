function [] = refreshBoundProps(obj)
% refreshBoundProps  Rebuild every bound propeller against the current
%                    strip grid and the current Tilt.
%
%   obj.refreshBoundProps()
%
%   discretise() rebuilds obj.Strips from scratch, which invalidates the
%   Mask / Cover / Swirl arrays cached in obj.BoundProps: they are sized to
%   the OLD strip count and were computed in the OLD tilted frame. This
%   replays each stored Spec through addBoundProps so the props follow the
%   surface.
%
%   Call this after ANY re-discretise - tiltSurf() already does.
%
%   No-op if there are no bound props. Definition order is preserved, so
%   BoundProps(j) refers to the same physical prop before and after.

if ~isprop(obj,'BoundProps') || isempty(obj.BoundProps)
    return;
end

if isempty(obj.Strips) || ~isfield(obj.Strips,'n')
    error('refreshBoundProps:NoStrips', ...
        'obj.Strips is empty - call discretise() before refreshing props.');
end

specs = {obj.BoundProps.Spec};      % 1 x nProps cell of arg lists

% clear first: addBoundProps appends, so replaying onto a populated
% list would duplicate every prop
obj.BoundProps = struct([]);
obj.Propwash   = [];
obj.PropCover  = [];

for j = 1:numel(specs)
    obj.addBoundProps(specs{j}{:});
end
end