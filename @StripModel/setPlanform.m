function obj = setPlanform(obj, opts)
% setPlanform  Set any subset of planform parameters by name-value.
%   wing.setPlanform('Span',0.5,'RootChord',0.18,'TaperRatio',0.7, ...
%                    'TwistTip',-3,'Incidence',2,'SweepQC',10)
    arguments
        obj
        opts.Span        (1,1) double {mustBePositive}
        opts.RootChord   (1,1) double {mustBePositive}
        opts.TaperRatio  (1,1) double {mustBePositive}
        opts.TwistTip    (1,1) double
        opts.Incidence   (1,1) double
        opts.SweepQC     (1,1) double
    end
    f = fieldnames(opts);
    for k = 1:numel(f), obj.(f{k}) = opts.(f{k}); end
end
