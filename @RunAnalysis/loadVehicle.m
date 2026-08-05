function plane = loadVehicle(obj)
% Runs the config script in this helper's workspace and hands back `plane`.
run(obj.ConfigScriptName);
if ~exist('plane', 'var')
    error("RunAnalysis:NoVehicle", ...
        "Config '%s' did not define a 'plane' struct.", obj.ConfigScriptName);
end

end