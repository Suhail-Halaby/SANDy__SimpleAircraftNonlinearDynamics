function S = discretise(obj, opts)
% discretise  Build spanwise strip geometry (pure geometry, no aero).
    arguments
        obj
        opts.nStrips (1,1) double {mustBeInteger,mustBePositive} = 20
        opts.Spacing (1,1) string {mustBeMember(opts.Spacing,["uniform","cosine"])} = "uniform"
    end
    n = opts.nStrips;

    % spanwise edges in eta = [0,1]
    switch opts.Spacing
        case "uniform"
            edges = linspace(0,1,n+1);
        case "cosine"
            edges = 0.5*(1 - cos(linspace(0,pi,n+1)));
    end
    eta  = 0.5*(edges(1:end-1) + edges(2:end));   % strip centres
    dEta = diff(edges);

    % orientation
    G  = deg2rad(obj.Dihedral);
    Lm = deg2rad(obj.SweepQC);
    spanHat = [0; obj.SpanSide*cos(G); -sin(G)];   % unit span vector, body axes

    % per-strip geometry
    s     = eta * obj.Span;                                   % span coord [m]
    ds    = dEta * obj.Span;                                  % strip width [m]
    chord = obj.RootChord * (1 - (1-obj.TaperRatio)*eta);     % local chord [m]
    twist = obj.Incidence + obj.Tilt + obj.TwistTip*eta;      % local incidence [deg]
    area  = chord .* ds;                                      % strip area [m^2]

    % quarter-chord position of each strip (UNTILTED), body axes.
    % Sweep applied along body -X (aft).
    qc = obj.RootQC + spanHat.*s + [-tan(Lm)*s; zeros(1,n); zeros(1,n)];

    % ---- per-strip local frame (UNTILTED) ----
    % chordHat: TE->LE (forward), rotated by incidence+twist about the span
    % axis. SpanSide keeps +incidence = LE-up on both halves.
    % NOTE: obj.Tilt is deliberately EXCLUDED here - it is applied once
    % below via Ry. Including it here would double-count the tilt.
    ang0      = obj.SpanSide * deg2rad(obj.Incidence + obj.TwistTip*eta);   % 1 x n
    chordHat  = [cos(ang0); spanHat(3)*sin(ang0); -spanHat(2)*sin(ang0)];   % 3 x n
    normalHat = cross(repmat(spanHat,1,n), chordHat);                       % 3 x n

    % ---- tilt: single rotation about body Y through RootQC+TiltOffset ----
    % rotate the QC line and the frame directions together => no doubling
    if obj.Tilt ~= 0
        t  = deg2rad(obj.Tilt);
        Ry = [cos(t) 0 sin(t); 0 1 0; -sin(t) 0 cos(t)];   % +Y, +t => LE up
        p  = obj.RootQC + obj.TiltOffset;
        qc        = p + Ry*(qc - p);
        chordHat  = Ry*chordHat;
        normalHat = Ry*normalHat;
        spanHat   = Ry*spanHat;     % [0;±1;0] for dihedral 0; exact in general
    end

    % map strips to control surfaces (membership by strip centre)
    ctrlIdx       = zeros(1,n);
    ctrlChordFrac = zeros(1,n);
    for i = 1:n
        for j = 1:numel(obj.Controls)
            c = obj.Controls(j);
            if eta(i) >= c.EtaInner && eta(i) <= c.EtaOuter
                ctrlIdx(i)       = j;
                ctrlChordFrac(i) = c.ChordFraction;
                break;   % first match wins
            end
        end
    end

    % reference quantities
    Sref = sum(area);
    MAC  = sum(chord.^2 .* ds) / Sref;   % mean aerodynamic chord
    mc = qc - 0.25 * chord .* chordHat;

    S = struct();
    S.n             = n;
    S.spacing       = opts.Spacing;
    S.edges         = edges;     % 1 x (n+1) eta edges
    S.eta           = eta;       % 1 x n      strip-centre span fraction
    S.s             = s;         % 1 x n      span coord [m]
    S.ds            = ds;        % 1 x n      strip width [m]
    S.chord         = chord;     % 1 x n      local chord [m]
    S.twist         = twist;     % 1 x n      local incidence [deg] (also in chordHat/normalHat)
    S.area          = area;      % 1 x n      strip area [m^2]
    S.qc            = qc;        % 3 x n      QC position, body axes [m]
    S.spanHat       = spanHat ./ vecnorm(spanHat);   % 3 x 1      span unit vector
    S.chordHat      = chordHat ./ vecnorm(chordHat);  % 3 x n      chord unit vector (TE->LE, forward)
    S.normalHat     = normalHat ./ vecnorm(normalHat); % 3 x n      surface normal (lift direction)
    S.ctrlIdx       = ctrlIdx;   % 1 x n      control index (0 = none)
    S.ctrlChordFrac = ctrlChordFrac;
    S.Sref          = Sref;      % panel reference area [m^2]
    S.MAC           = MAC;       % mean aerodynamic chord [m]
    S.b             = obj.Span;  % panel span [m]
    S.Colloc        = qc;        % 3 x n  quarter-chord = collocation point [m]
    S.mc  =  mc;                 % 3 x n  mid-chord [m]

    obj.Strips = S;
end
