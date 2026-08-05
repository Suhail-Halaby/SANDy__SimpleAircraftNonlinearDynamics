function fig = StripPlotter(varargin)
% STRIPPLOTTER  Wireframe plotter for one or more StripModel lifting surfaces.
%
%   StripPlotter(wingR)                       % single surface
%   StripPlotter(wingR, wingL, fin)           % several surfaces, shared axes
%   StripPlotter(surfaceArray)                % a StripModel array
%   StripPlotter(..., 'Strips',false, 'QC',true, 'Title',"my aircraft")
%
%   fig = StripPlotter(...)  returns the figure handle.
%
% Draws, per surface, in body axes (X fwd, Y stbd, Z down) relative to CG:
%   - planform outline (LE / TE / root & tip chords)   [always]
%   - spanwise strip boundaries                        [toggle: key 's']
%   - control-surface regions + hinge lines            [toggle: key 'c']
%   - quarter-chord line                               [toggle: key 'q']
%   - bound propeller disks + hub dots                 [toggle: key 'p']
%
% Geometry is reconstructed from the object's planform properties, so it works
% whether or not discretise() has been called. If it has, the surface's own
% strip edges (and spacing) are used; otherwise a default grid is shown.
%
% Propeller disks are read from obj.BoundProps (populated by addBoundProps).
% Unlike the planform, they are NOT reconstructed - a surface with no bound
% props simply draws none.
%
% Name-value options:
%   'Strips'   (true)   initial visibility of strip boundaries
%   'Controls' (true)   initial visibility of control surfaces
%   'QC'       (false)  initial visibility of quarter-chord line
%   'Props'    (true)   initial visibility of propeller disks
%   'Title'    ("Strip model")
%   'Triad'    (true)   draw body-axis triad at CG
%
% ---- Debug notes ----
% Each group of graphics objects is Tagged ('strips'/'controls'/'qc'/'props'),
% so the toggle state lives on the objects, not in a separate bookkeeping
% struct. To poke at a group from the command line:
%   findobj(gca,'Tag','props')

    [models, o] = parseInputs(varargin);

    fig = figure('Color','w','Name','StripPlotter');
    ax  = axes('Parent',fig); hold(ax,'on');
    C   = lines(numel(models));

    legH = gobjects(1,numel(models));
    for i = 1:numel(models)
        legH(i) = drawSurface(ax, models(i), C(i,:));
    end

    drawCGandTriad(ax, o.Triad);
    styleAxes(ax, o.Title);
    legend(legH, 'Location','northeastoutside');

    % initial visibility (state lives on the tagged objects)
    setGroup(ax, 'strips',   o.Strips);
    setGroup(ax, 'controls', o.Controls);
    setGroup(ax, 'qc',       o.QC);
    setGroup(ax, 'props',    o.Props);

    set(fig, 'KeyPressFcn', @(~,evt) onKey(ax, evt));
end

% ===================================================================
% Input parsing
% ===================================================================
function [models, o] = parseInputs(args)
    models = StripModel.empty;
    k = 1;
    while k <= numel(args) && isa(args{k},'StripModel')
        models = [models, reshape(args{k},1,[])]; %#ok<AGROW>
        k = k + 1;
    end
    if isempty(models)
        error('StripPlotter:noModels','Pass at least one StripModel object.');
    end

    p = inputParser;
    p.addParameter('Strips',   true,  @(x)islogical(x)||isnumeric(x));
    p.addParameter('Controls', true,  @(x)islogical(x)||isnumeric(x));
    p.addParameter('QC',       false, @(x)islogical(x)||isnumeric(x));
    p.addParameter('Props',    true,  @(x)islogical(x)||isnumeric(x));
    p.addParameter('Title',    "Strip model", @(x)isstring(x)||ischar(x));
    p.addParameter('Triad',    true,  @(x)islogical(x)||isnumeric(x));
    p.parse(args{k:end});
    o = p.Results;
end

% ===================================================================
% Drawing  (one surface = planform + strips + QC + controls + props)
% ===================================================================
function legHandle = drawSurface(ax, m, col)
    legHandle = drawPlanform(ax, m, col);   % returns the handle used for the legend
    drawStrips(ax, m, col);
    drawQuarterChord(ax, m, col);
    drawControls(ax, m, col);
    drawProps(ax, m, col);
end

function h = drawPlanform(ax, m, col)
    e  = linspace(0,1,40);
    LE = surfPoint(m, e,     zeros(size(e)));   % xc = 0  (leading edge)
    TE = surfPoint(m, e,     ones(size(e)));    % xc = 1  (trailing edge)
    rc = surfPoint(m, [0 0], [0 1]);            % root chord
    tc = surfPoint(m, [1 1], [0 1]);            % tip chord

    h = plot3(ax, LE(1,:),LE(2,:),LE(3,:), '-','Color',col,'LineWidth',2, ...
              'DisplayName',m.Name);
    plot3(ax, TE(1,:),TE(2,:),TE(3,:), '-','Color',col,'LineWidth',2,'HandleVisibility','off');
    plot3(ax, rc(1,:),rc(2,:),rc(3,:), '-','Color',col,'LineWidth',2,'HandleVisibility','off');
    plot3(ax, tc(1,:),tc(2,:),tc(3,:), '-','Color',col,'LineWidth',2,'HandleVisibility','off');
end

function drawStrips(ax, m, col)
    if isfield(m.Strips,'edges'), edges = m.Strips.edges;
    else,                         edges = linspace(0,1,11);   % default grid
    end
    for ee = edges
        P = surfPoint(m, [ee ee], [0 1]);       % chordwise line at this eta
        plot3(ax, P(1,:),P(2,:),P(3,:), '-','Color',[col 0.5],'LineWidth',0.5, ...
              'Tag','strips','HandleVisibility','off');
    end
end

function drawQuarterChord(ax, m, col)
    QC = surfPoint(m, linspace(0,1,40), 0.25*ones(1,40));
    plot3(ax, QC(1,:),QC(2,:),QC(3,:), '--','Color',col*0.6,'LineWidth',1, ...
          'Tag','qc','HandleVisibility','off');
end

function drawControls(ax, m, col)
    for j = 1:numel(m.Controls)
        c   = m.Controls(j);
        eC  = linspace(c.EtaInner, c.EtaOuter, 20);
        xch = 1 - c.ChordFraction;                  % hinge chord fraction
        Hng = surfPoint(m, eC, xch*ones(size(eC)));  % hinge line
        Ted = surfPoint(m, eC, ones(size(eC)));      % trailing edge
        poly = [Hng, fliplr(Ted)];                   % closed region

        patch(ax, 'XData',poly(1,:),'YData',poly(2,:),'ZData',poly(3,:), ...
              'FaceColor',col,'FaceAlpha',0.30,'EdgeColor','none', ...
              'Tag','controls','HandleVisibility','off');
        plot3(ax, Hng(1,:),Hng(2,:),Hng(3,:), '-','Color',col,'LineWidth',1.5, ...
              'Tag','controls','HandleVisibility','off');
    end
end

% ===================================================================
% Bound propellers: disk in the plane of rotation + hub dot
% ===================================================================
function drawProps(ax, m, col)
% Disks come straight from obj.BoundProps (Hub, Axis, Radius) - no
% reconstruction, so this is a silent no-op if addBoundProps was never called.
    if ~isprop(m,'BoundProps') || isempty(m.BoundProps), return; end

    th   = linspace(0, 2*pi, 60);
    dark = col*0.5;

    for j = 1:numel(m.BoundProps)
        P = m.BoundProps(j);

        a = P.Axis(:);
        if norm(a) == 0, continue; end
        a = a / norm(a);
        [u,v] = orthoBasis(a);                       % disk plane, normal = a

        ring = P.Hub(:) + P.Radius*(u*cos(th) + v*sin(th));   % 3 x 60

        patch(ax, 'XData',ring(1,:), 'YData',ring(2,:), 'ZData',ring(3,:), ...
              'FaceColor',col, 'FaceAlpha',0.15, 'EdgeColor',dark, ...
              'LineWidth',1.2, 'Tag','props', 'HandleVisibility','off');

        plot3(ax, P.Hub(1), P.Hub(2), P.Hub(3), '.', 'Color',dark, ...
              'MarkerSize',18, 'Tag','props', 'HandleVisibility','off');
    end
end

function [u,v] = orthoBasis(a)
% Two unit vectors spanning the plane normal to unit vector a.
% Seed is chosen away from a to keep the cross product well conditioned.
    if abs(a(1)) < 0.9, t = [1;0;0]; else, t = [0;1;0]; end
    u = cross(a,t);  u = u / norm(u);
    v = cross(a,u);                     % already unit: a,u orthonormal
end

function drawCGandTriad(ax, showTriad)
    plot3(ax, 0,0,0, 'k+','MarkerSize',10,'LineWidth',1.5,'HandleVisibility','off');
    if showTriad
        L = 0.18 * max([diff(xlim(ax)) diff(ylim(ax)) diff(zlim(ax)) eps]);
        triad(ax, L);
    end
end

% ===================================================================
% Axes cosmetics / view  (single place to tweak orientation)
% ===================================================================
function styleAxes(ax, ttl)
    axis(ax,'equal'); grid(ax,'on'); box(ax,'on');
    set(ax,'ZDir','reverse');                 % body +Z is down -> plot downward
    set(gca,'YDir','reverse')
    xlabel(ax,'x  (fwd) [m]'); ylabel(ax,'y  (stbd) [m]'); zlabel(ax,'z  (down) [m]');
    view(ax,-135,20);                         % <-- your fixed azimuth goes here
    enableDefaultInteractivity(ax);
    title(ax, ttl);
    subtitle(ax, 'keys:  s = strips   c = controls   q = quarter-chord   p = props');
end

% ===================================================================
% Visibility toggles  (state read straight off the tagged objects)
% ===================================================================
function onKey(ax, evt)
    switch evt.Key
        case 's', toggleGroup(ax, 'strips');
        case 'c', toggleGroup(ax, 'controls');
        case 'q', toggleGroup(ax, 'qc');
        case 'p', toggleGroup(ax, 'props');
    end
end

function toggleGroup(ax, tag)
    h = findobj(ax, 'Tag', tag);
    if isempty(h), return; end
    isOn = strcmp(get(h(1),'Visible'), 'on');
    set(h, 'Visible', ~isOn);
end

function setGroup(ax, tag, tf)
    set(findobj(ax,'Tag',tag), 'Visible', logical(tf));   % set([],...) is a no-op
end

% ===================================================================
% Geometry: 3-D body-axis point on a surface
% ===================================================================
function P = surfPoint(m, eta, xc)
% 3-D body-axis point on surface m at span fraction eta (0..1, root->tip)
% and chord fraction xc (0 = LE, 0.25 = QC, 1 = TE). eta, xc are 1xK.
    G  = deg2rad(m.Dihedral);
    Lm = deg2rad(m.SweepQC);
    spanHat = [0; m.SpanSide*cos(G); -sin(G)];          % root->tip unit vector

    s     = eta * m.Span;                                % span coord [m]
    chord = m.RootChord .* (1 - (1-m.TaperRatio).*eta);  % local chord [m]
    inc   = deg2rad(m.Incidence + m.Tilt + m.TwistTip.*eta);  % local incidence [rad]

    % quarter-chord position (sweep applied aft along body -X)
    qc = m.RootQC + spanHat.*s + [-tan(Lm).*s; zeros(1,numel(eta)); zeros(1,numel(eta))];

    % tilt: rotate the QC line about body Y through RootQC+TiltOffset
    % (chord tilt is handled via 'inc' above, so rotate QC only - matches discretise)
    if m.Tilt ~= 0
        t  = deg2rad(m.Tilt);
        Ry = [cos(t) 0 sin(t); 0 1 0; -sin(t) 0 cos(t)];
        p  = m.RootQC + m.TiltOffset;
        qc = p + Ry*(qc - p);
    end

    % chordwise unit vector (TE->LE forward), incidence about the span axis.
    % SpanSide keeps positive incidence = LE-up on both halves.
    ang = m.SpanSide .* inc;
    chordHat = [cos(ang); spanHat(3)*sin(ang); -spanHat(2)*sin(ang)];

    % offset from QC: xc=0.25 gives QC, xc=0 gives +0.25c fwd, xc=1 gives -0.75c aft
    P = qc + ((0.25 - xc).*chord) .* chordHat;
end

% ===================================================================
% Body-axis triad at the CG
% ===================================================================
function triad(ax, L)
    axisVec = [L 0 0; 0 L 0; 0 0 L];
    cols    = [0.85 0 0; 0 0.7 0; 0 0 0.85];
    labs    = {'x','y','z'};
    for a = 1:3
        v = axisVec(a,:);
        quiver3(ax, 0,0,0, v(1),v(2),v(3), 0, 'Color',cols(a,:), ...
                'LineWidth',1.5, 'MaxHeadSize',0.6, 'HandleVisibility','off');
        text(ax, 1.15*v(1), 1.15*v(2), 1.15*v(3), labs{a}, ...
             'Color',cols(a,:), 'FontWeight','bold');
    end
end