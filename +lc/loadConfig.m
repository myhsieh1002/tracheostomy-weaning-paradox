function cfg = loadConfig(configPath)
%LOADCONFIG Read the Model 2 configuration file into a struct.
%
%   cfg = lc.loadConfig() reads config/params_m2.json.
%
%   As in Model 1, {value, source} wrappers are flattened to plain scalars
%   in cfg.dynamics while provenance is preserved in cfg.sources and the
%   sensitivity ranges in cfg.ranges.
%
%   See also wob.loadConfig, lc.params

arguments
    configPath (1,:) char = fullfile(wob.projectRoot, 'config', 'params_m2.json')
end

if ~isfile(configPath)
    error('lc:loadConfig:fileNotFound', 'Config file not found: %s', configPath);
end

raw = jsondecode(fileread(configPath));

cfg = raw;
cfg.configPath = configPath;
cfg.sources = struct();
cfg.ranges  = struct();

names = fieldnames(raw.dynamics);
for k = 1:numel(names)
    name = names{k};
    entry = raw.dynamics.(name);
    if isstruct(entry) && isfield(entry, 'value')
        cfg.dynamics.(name) = entry.value;
        cfg.sources.(name)  = entry.source;
        if isfield(entry, 'range')
            cfg.ranges.(name) = reshape(entry.range, 1, []);
        end
    end
end

cfg.capacity_grid.C_max = reshape(raw.capacity_grid.C_max, 1, []);
cfg.continuation.L_range = reshape(raw.continuation.L_range, 1, []);
end
