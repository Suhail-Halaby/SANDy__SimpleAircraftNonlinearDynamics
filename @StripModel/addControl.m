function obj = addControl(obj, opts)
% addControl  Append a trailing-edge control surface.
%   EtaInner/EtaOuter : span fraction [0..1] it occupies
%   ChordFraction     : control chord / local chord [0..1]
%   wing.addControl('Name',"elevon",'EtaInner',0.4,'EtaOuter',0.95, ...
%                   'ChordFraction',0.25,'Limits',[-25 25])
    arguments
        obj
        opts.Name          (1,1) string = "ctrl" + string(numel(obj.Controls)+1)
        opts.EtaInner      (1,1) double {mustBeInRange(opts.EtaInner,0,1)}
        opts.EtaOuter      (1,1) double {mustBeInRange(opts.EtaOuter,0,1)}
        opts.ChordFraction (1,1) double {mustBeInRange(opts.ChordFraction,0,1)}
        opts.Limits        (1,2) double = deg2rad([-30 30])  % deflection limits [rad]
        opts.Deflection    (1,1) double = 0          % current command [rad]
        opts.Direction (1,1) double = 1;
    end

    if opts.EtaInner >= opts.EtaOuter
        error('StripModel:badControl','EtaInner must be < EtaOuter.');
    end

    c.Name          = opts.Name;
    c.EtaInner      = opts.EtaInner;
    c.EtaOuter      = opts.EtaOuter;
    c.ChordFraction = opts.ChordFraction;
    c.Deflection    = opts.Deflection;
    c.Limits        = opts.Limits;
    c.Direction     = opts.Direction; 
    obj.Controls(end+1) = c;
end
