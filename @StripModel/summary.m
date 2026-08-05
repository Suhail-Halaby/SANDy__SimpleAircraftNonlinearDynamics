function summary(obj)
% summary  Print a short text description of the panel to the command window.
    fprintf('StripModel "%s"\n', obj.Name);
    fprintf('  span=%.3f m  rootC=%.3f m  taper=%.3f  tipC=%.3f m\n', ...
        obj.Span, obj.RootChord, obj.TaperRatio, obj.RootChord*obj.TaperRatio);
    fprintf('  incidence=%.2f  twistTip=%.2f  sweepQC=%.2f  dihedral=%.2f deg  side=%+d\n', ...
        obj.Incidence, obj.TwistTip, obj.SweepQC, obj.Dihedral, obj.SpanSide);
    fprintf('  rootQC=[% .3f % .3f % .3f] m (body, rel CG)\n', obj.RootQC);
    fprintf('  controls: %d\n', numel(obj.Controls));
    for j = 1:numel(obj.Controls)
        c = obj.Controls(j);
        fprintf('    [%d] %-10s eta[%.2f-%.2f]  cf=%.2f  lim[%g %g]\n', ...
            j, c.Name, c.EtaInner, c.EtaOuter, c.ChordFraction, c.Limits(1), c.Limits(2));
    end
    if isfield(obj.Strips,'Sref')
        fprintf('  discretised: %d strips  S=%.4f m^2  MAC=%.3f m\n', ...
            obj.Strips.n, obj.Strips.Sref, obj.Strips.MAC);
    end
end
