function setupPaths()
root = fileparts(mfilename('fullpath'));
addpath(root);
addpath(genpath(fullfile(root,'Utilities')));
addpath(genpath(fullfile(root,'Examples')));
end