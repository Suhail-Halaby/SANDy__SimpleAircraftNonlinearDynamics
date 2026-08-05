function [] = addBoundProps(obj, opts)
% addBoundProps  Attach a wing-mounted (bound) propeller and flag the strips
%                that lie in its immediate wake.
%
%   obj.addBoundProps(Radius=..., Eta=..., ChordFrac=...)
%
%   Appends a propeller definition to obj.BoundProps and updates
%   obj.propwash - a 1 x n array of ones/zeros, one entry per strip, where
%   1 means the strip sits inside the immediate slipstream of at least one
%   bound propeller.
%
%   Requires discretise() to have been called first (uses obj.Strips).
%
%   Name-value options
%     Radius     prop radius [m]                                (required)
%     Eta        hub location along the wing, eta in [0,1]      (required)
%     ChordFrac  disk-plane position as x/c from the LE.
%                0 = LE, 1 = TE, negative = ahead of the LE.
%                Default -0.25 (tractor, quarter-chord ahead of LE).
%     Axis       thrust axis, UNTILTED body axes (see FRAMES below).
%                Default [0;0;0] => use the local chordHat at the hub
%                station, i.e. aligned with the twisted chord line.
%     Offset     extra hub offset from the chord plane, UNTILTED body axes [m].
%     AxisFixed  false (default): the prop is bolted to this surface, so Axis
%                and Offset rotate with obj.Tilt.
%                true: Axis and Offset are taken as final body-axis vectors
%                and are NOT rotated by Tilt (body/fuselage-mounted prop).
%     Spin       +1 / -1, rotation sense about Axis (right-hand about +Axis).
%                Only used to set the per-strip swirl sign.
%     Name       optional label.
%     nSample    spanwise sample points per strip used for coverage. Default 9.
%     MinCover   strip flagged when coverage > MinCover. Default 0
%                (any overlap at all counts).
%
%   ---------------------------------------------------------------------
%   FRAMES
%   Mirrors discretise(): tilt is applied ONCE, via Ry about body Y through
%   RootQC+TiltOffset. The default Axis (chordHat) and the hub station (qc)
%   already carry the tilt because they come from obj.Strips, which
%   discretise() has already rotated - so Ry is applied ONLY to an
%   explicitly supplied Axis/Offset. Do not pre-tilt those yourself.
%   ---------------------------------------------------------------------

%
%   Every call stores its own arguments in BoundProps(j).Spec so the prop can
%   be rebuilt against a new strip grid - see refreshBoundProps(), which
%   tiltSurf() and any re-discretise must call.
%
%   Wake model: the immediate wake is the semi-infinite cylinder of radius
%   Radius, aligned with Axis, extending downstream (-Axis) from the disk
%   plane. Contraction, skew and slipstream decay are deliberately NOT
%   modelled here - this is a pure geometric membership test, in the same
%   spirit as discretise().

    arguments
        obj
        opts.Radius    (1,1) double {mustBePositive}
        opts.Eta       (1,1) double {mustBeGreaterThanOrEqual(opts.Eta,0), ...
                                     mustBeLessThanOrEqual(opts.Eta,1)}
        opts.ChordFrac (1,1) double = -0.25
        opts.Axis      (3,1) double = [0;0;0]
        opts.Offset    (3,1) double = [0;0;0]
        opts.AxisFixed (1,1) logical = false
        opts.Spin      (1,1) double {mustBeMember(opts.Spin,[-1 1])} = 1
        opts.Name      (1,1) string = ""
        opts.nSample   (1,1) double {mustBeInteger,mustBePositive} = 9
        opts.MinCover  (1,1) double {mustBeGreaterThanOrEqual(opts.MinCover,0), ...
                                     mustBeLessThanOrEqual(opts.MinCover,1)} = 0
    end

    if isempty(obj.Strips) || ~isfield(obj.Strips,'n')
        error('addBoundProps:NoStrips', ...
              'Call discretise() before addBoundProps().');
    end
    S = obj.Strips;
    n = S.n;

    % ---- verbatim spec, for replay against a rebuilt strip grid -----------
    spec = {'Radius',    opts.Radius, ...
            'Eta',       opts.Eta, ...
            'ChordFrac', opts.ChordFrac, ...
            'Axis',      opts.Axis, ...
            'Offset',    opts.Offset, ...
            'AxisFixed', opts.AxisFixed, ...
            'Spin',      opts.Spin, ...
            'Name',      opts.Name, ...
            'nSample',   opts.nSample, ...
            'MinCover',  opts.MinCover};

    % ---- hub station geometry -------------------------------------------
    % qc is linear in eta (rigid rotation of a linear span sweep), so
    % linear interp/extrap is exact. chordHat is not linear when TwistTip
    % ~= 0, so it is interpolated then re-normalised.
    qcHub    = stationVec(S, S.qc,       opts.Eta);
    chordHub = stationVec(S, S.chordHat, opts.Eta);
    chordHub = chordHub / norm(chordHub);
    spanHat  = S.spanHat / norm(S.spanHat);

    cHub = obj.RootChord * (1 - (1-obj.TaperRatio)*opts.Eta);   % local chord

    % ---- tilt for explicitly supplied vectors only ------------------------
    axisIn   = opts.Axis;
    offsetIn = opts.Offset;
    if obj.Tilt ~= 0 && ~opts.AxisFixed
        t  = deg2rad(obj.Tilt);
        Ry = [cos(t) 0 sin(t); 0 1 0; -sin(t) 0 cos(t)];
        if any(axisIn ~= 0), axisIn = Ry*axisIn; end   % zeros => chordHat, already tilted
        offsetIn = Ry*offsetIn;
    end

    if all(axisIn == 0)
        axisHat = chordHub;                 % thrust along the local chord line
    else
        axisHat = axisIn / norm(axisIn);
    end

    % chordHat points TE -> LE, and qc sits 25% aft of the LE, so a point at
    % fraction f of the chord from the LE is qc + (0.25 - f)*c*chordHat.
    hub = qcHub + (0.25 - opts.ChordFrac)*cHub*chordHub + offsetIn;

    % ---- per-strip coverage ---------------------------------------------
    dEta  = diff(S.edges);
    u     = linspace(-0.5, 0.5, opts.nSample);
    cover = zeros(1,n);
    swirl = zeros(1,n);

    for i = 1:n
        etaS = S.eta(i) + u*dEta(i);        % sample across the strip width
        Pw   = stationVec(S, S.qc, etaS);   % 3 x nSample, collocation line
        d    = Pw - hub;

        axial  = axisHat.' * d;                          % 1 x k, +ve = ahead of disk
        radial = vecnorm(d - axisHat*axial, 2, 1);       % 1 x k

        inside   = (axial <= 0) & (radial <= opts.Radius);
        cover(i) = mean(inside);

        % swirl sign: up-going vs down-going blade side of the hub
        yOff     = spanHat.' * (S.Colloc(:,i) - hub);
        swirl(i) = -opts.Spin * sign(yOff);
    end

    mask  = double(cover > opts.MinCover);
    swirl = swirl .* mask;

    % ---- store -----------------------------------------------------------
    % Built field-by-field: struct('Spec',spec) with a cell value would
    % expand into a struct ARRAY, which is not what we want.
    Pdef           = struct();
    Pdef.Name      = opts.Name;
    Pdef.Radius    = opts.Radius;
    Pdef.Eta       = opts.Eta;
    Pdef.ChordFrac = opts.ChordFrac;
    Pdef.Hub       = hub;         % 3 x 1 body axes [m], tilt applied
    Pdef.Axis      = axisHat;     % 3 x 1 unit, points forward, tilt applied
    Pdef.Spin      = opts.Spin;
    Pdef.Mask      = mask;        % 1 x n ones/zeros
    Pdef.Cover     = cover;       % 1 x n fraction of strip washed
    Pdef.Swirl     = swirl;       % 1 x n  +1/-1/0
    Pdef.Spec      = spec;        % 1 x 2m cell, replay arguments

    if isempty(obj.BoundProps)
        obj.BoundProps = Pdef;
    else
        obj.BoundProps(end+1) = Pdef;
    end

    % combined across every bound prop on this surface
    allMask       = vertcat(obj.BoundProps.Mask);        % nProps x n
    obj.Propwash  = double(any(allMask, 1));             % 1 x n ones/zeros
    obj.PropCover = min(1, sum(vertcat(obj.BoundProps.Cover), 1));
end

% =========================================================================
function V = stationVec(S, A, etaQ)
% Interpolate a 3 x n strip quantity onto arbitrary eta stations.
    if S.n >= 2
        V = interp1(S.eta(:), A.', etaQ(:), 'linear', 'extrap').';
    else
        V = repmat(A(:,1), 1, numel(etaQ));
    end
end