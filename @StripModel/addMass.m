function [] = addMass(obj, opts)
% addMass  Build the inertial contribution of this surface: total mass, CG,
%          and inertia tensor about that CG, in body axes.
%
%   obj.addMass(AreaDensity=2.1)
%   obj.addMass(Shape="flatplate", AreaDensity=2.1, ...
%               AttachedPointMass=[0.25 0.25], PointMassLoc=[x1 x2; y1 y2; z1 z2])
%
%   Requires discretise() to have been called first (uses obj.Strips).
%
%   Name-value options
%     Shape              "lumpmass"  : each strip is a point mass at its
%                                      mid-chord (no self-inertia).
%                        "flatplate" : each strip additionally carries the
%                                      inertia of a thin rectangular plate of
%                                      chord x ds, in its own local frame.
%     AreaDensity        surface mass per unit area [kg/m^2]. Default 0.
%     AttachedPointMass  1 x p masses [kg].
%     PointMassLoc       3 x p positions, UNTILTED body axes [m] (see FRAMES).
%     LocFixed           false (default): point masses are bolted to this
%                                      surface and rotate with obj.Tilt.
%                        true:  PointMassLoc is taken as final body axes and
%                               does NOT follow the tilt.
%
%   Sets obj.PartMass, obj.PartCG, obj.Inertias, and stores the call
%   arguments in obj.MassSpec so the mass model can be rebuilt against a new
%   strip grid - see refreshMass(), which tiltSurf() calls.
%
%   ---------------------------------------------------------------------
%   FRAMES
%   Mirrors discretise(): tilt is applied ONCE, via Ry about body Y through
%   RootQC+TiltOffset. Strip quantities (mc, chordHat, normalHat) come from
%   obj.Strips and already carry the tilt. Ry is therefore applied ONLY to
%   PointMassLoc - do not pre-tilt those yourself.
%
%   NOTE: PointMassLoc holds POSITIONS, so it is rotated about the hinge
%   point p = RootQC + TiltOffset as p + Ry*(loc - p), not by a bare Ry.
%   ---------------------------------------------------------------------
%   Add to your classdef properties block:
%   ---------------------------------------------------------------------

    arguments
        obj
        opts.Shape (1,1) string {mustBeMember(opts.Shape,["lumpmass" "flatplate"])} = "flatplate"
        opts.AreaDensity (1,1) double {mustBeNonnegative} = 0
        opts.AttachedPointMass (1,:) double {mustBeNonnegative} = []
        opts.PointMassLoc (3,:) double = []
        opts.LocFixed (1,1) logical = false
    end

    if numel(opts.AttachedPointMass) ~= size(opts.PointMassLoc,2)
        error('StripModel:addMass','Each point mass needs one xyz column in PointMassLoc.');
    end
    if isempty(obj.Strips) || ~isfield(obj.Strips,'n')
        error('StripModel:addMass','Call discretise() before addMass().');
    end

    % ---- verbatim spec, for replay against a rebuilt strip grid -----------
    obj.MassSpec = {'Shape',             opts.Shape, ...
                    'AreaDensity',       opts.AreaDensity, ...
                    'AttachedPointMass', opts.AttachedPointMass, ...
                    'PointMassLoc',      opts.PointMassLoc, ...
                    'LocFixed',          opts.LocFixed};

    S = obj.Strips;
    n = S.n;

    % ---- distributed mass: one lump per strip, at its mid-chord -----------
    % (mc, not the collocation point - the strip's mass centroid is at
    %  mid-chord, the collocation point is at the quarter-chord.)
    stripMass = S.area * opts.AreaDensity;   % 1 x n  [kg]
    loc       = S.mc;                        % 3 x n  body axes [m], tilted

    % ---- point masses: rotate about the tilt hinge unless held fixed ------
    pmLoc = opts.PointMassLoc;
    if obj.Tilt ~= 0 && ~opts.LocFixed && ~isempty(pmLoc)
        t  = deg2rad(obj.Tilt);
        Ry = [cos(t) 0 sin(t); 0 1 0; -sin(t) 0 cos(t)];
        p  = obj.RootQC + obj.TiltOffset;
        pmLoc = p + Ry*(pmLoc - p);          % POSITION: rotate about p
    end

    massAll = [stripMass, opts.AttachedPointMass];   % 1 x (n+p)
    locAll  = [loc,       pmLoc];                    % 3 x (n+p)

    mTotal = sum(massAll);
    if mTotal <= 0
        error('StripModel:addMass', ...
              ['Total mass is zero - CG and inertia are undefined. ' ...
               'Set AreaDensity > 0 and/or supply AttachedPointMass.']);
    end

    partCg = (locAll * massAll.') / mTotal;   % 3 x 1  mass-weighted centroid

    % ---- inertia tensor about the CG -------------------------------------
    % Point-mass (parallel-axis) sum over every lump, for both shapes.
    I = zeros(3);
    for k = 1:numel(massAll)
        r = locAll(:,k) - partCg;                       % offset from CG
        I = I + massAll(k) * ((r.'*r)*eye(3) - r*r.');
    end

    if opts.Shape == "flatplate"
        % Each strip also carries its own thin-plate inertia about its
        % centroid. Local frame: e1 = chord (extent c), e2 = span (extent
        % ds), e3 = normal (zero thickness).
        %   I11 = m*ds^2/12,  I22 = m*c^2/12,  I33 = m*(c^2+ds^2)/12
        % Columns are re-normalised here because discretise() divides the
        % 3 x n chordHat/normalHat by norm() of the whole matrix, which is
        % a matrix 2-norm, not a per-column norm.
        e2 = S.spanHat / norm(S.spanHat);
        for i = 1:n
            m = stripMass(i);
            if m <= 0, continue; end
            c  = S.chord(i);
            ds = S.ds(i);

            e1 = S.chordHat(:,i);  e1 = e1 / norm(e1);
            e3 = S.normalHat(:,i); e3 = e3 / norm(e3);
            R  = [e1, e2, e3];                          % local -> body

            Iloc = diag([ m*ds^2/12, m*c^2/12, m*(c^2 + ds^2)/12 ]);
            I    = I + R*Iloc*R.';
        end
    end

    % ---- store on the object ---------------------------------------------
    obj.PartMass = mTotal;
    obj.PartCG   = partCg;
    obj.Inertias = I;
end