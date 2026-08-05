function m = mirror(obj, opts)
% mirror  Return a new StripModel reflected across the body XZ plane.
%
%   m = wing.mirror()              % counter-rotating prop
%   m = wing.mirror(RevProp=-1)    % co-rotating prop
%
% Specs are mirrored and replayed, not the derived fields: refreshBoundProps
% and refreshMass rebuild from the spec, so hand-edited Hub/Spin/Inertias
% would be overwritten on the first tiltSurf.

    arguments
        obj
        opts.RevProp (1,1) double {mustBeMember(opts.RevProp,[-1,1])} = 1
    end

    F = diag([1 -1 1]);
    m = StripModel(obj.Name + "_mirror");

    m.Span       = obj.Span;
    m.RootChord  = obj.RootChord;
    m.TaperRatio = obj.TaperRatio;
    m.TwistTip   = obj.TwistTip;
    m.Incidence  = obj.Incidence;
    m.SweepQC    = obj.SweepQC;
    m.Dihedral   = obj.Dihedral;
    m.Controls   = obj.Controls;

    m.SpanSide   = -obj.SpanSide;
    m.RootQC     = F * obj.RootQC;
    m.TiltOffset = F * obj.TiltOffset;
    m.Tilt       = obj.Tilt;

    m.MassSpec = obj.MassSpec;
    if ~isempty(m.MassSpec)
        m.MassSpec = specSet(m.MassSpec, 'PointMassLoc', ...
                             F * specGet(m.MassSpec,'PointMassLoc'));
    end

    propSpecs = cell(1, numel(obj.BoundProps));

    for j = 1:numel(obj.BoundProps)
        s = obj.BoundProps(j).Spec;
        s = specSet(s, 'Spin',   -opts.RevProp * specGet(s,'Spin'));
        s = specSet(s, 'Axis',    F * specGet(s,'Axis'));
        s = specSet(s, 'Offset',  F * specGet(s,'Offset'));
        s = specSet(s, 'Name',    specGet(s,'Name') + "_mirror");
        propSpecs{j} = s;
    end

    for k = 1:numel(obj.Controls)
        m.Controls.Direction = -obj.Controls(k).Direction;
    end


    if isfield(obj.Strips,'n')
        m.discretise('nStrips', obj.Strips.n, 'Spacing', obj.Strips.spacing);
        for j = 1:numel(propSpecs)
            m.addBoundProps(propSpecs{j}{:});
        end
        m.refreshMass();
    end
end

% ---- spec helpers: specs are {'Name',value,...} cells ----
function v = specGet(spec, name)
    k = find(strcmp(spec(1:2:end), name), 1);
    if isempty(k), v = []; else, v = spec{2*k}; end
end

function spec = specSet(spec, name, val)
    k = find(strcmp(spec(1:2:end), name), 1);
    if isempty(k), spec = [spec, {name, val}];
    else,          spec{2*k} = val;
    end
end