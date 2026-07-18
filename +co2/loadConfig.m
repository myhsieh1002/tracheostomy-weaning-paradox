function cfg = loadConfig(configPath)
%LOADCONFIG Read the Model 1b configuration into a struct.
%
%   cfg = co2.loadConfig() reads config/params_co2.json.
%
%   As in Models 1 and 2, {value, source} wrappers are flattened to plain
%   scalars while provenance is preserved in cfg.sources.
%
%   See also wob.loadConfig, lc.loadConfig

arguments
    configPath (1,:) char = fullfile(wob.projectRoot, 'config', 'params_co2.json')
end

if ~isfile(configPath)
    error('co2:loadConfig:fileNotFound', 'Config file not found: %s', configPath);
end

raw = jsondecode(fileread(configPath));
cfg = raw;
cfg.configPath = configPath;
cfg.sources = struct();

for section = ["metabolism", "deadspace", "stores", "targets"]
    if ~isfield(raw, section), continue; end
    names = fieldnames(raw.(section));
    for k = 1:numel(names)
        name = names{k};
        entry = raw.(section).(name);
        if isstruct(entry) && isfield(entry, 'value')
            cfg.(section).(name) = entry.value;
            cfg.sources.(name)   = entry.source;
        end
    end
end

cfg.metabolism.VCO2_grid  = reshape(raw.metabolism.VCO2_grid, 1, []);
cfg.deadspace.VD_VT_grid  = reshape(raw.deadspace.VD_VT_grid, 1, []);
end
