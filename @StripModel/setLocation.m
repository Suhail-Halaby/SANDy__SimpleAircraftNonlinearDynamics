function obj = setLocation(obj, opts)
% setLocation  Set root quarter-chord position and orientation.
%   wing.setLocation('RootQC',[0.05; 0.06; 0],'Dihedral',3,'SpanSide',+1)
    arguments
        obj
        opts.RootQC   (3,1) double
        opts.Dihedral (1,1) double
        opts.SpanSide (1,1) double {mustBeMember(opts.SpanSide,[-1 1])}
    end
    f = fieldnames(opts);
    for k = 1:numel(f), obj.(f{k}) = opts.(f{k}); end
end
