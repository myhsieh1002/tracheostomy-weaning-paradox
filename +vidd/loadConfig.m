function cfg = loadConfig(configPath)
%LOADCONFIG Read the Model 2b configuration into a struct.
%
%   cfg = vidd.loadConfig() reads config/params_vidd.json.
%
%   See also wob.loadConfig, lc.loadConfig, co2.loadConfig

arguments
    configPath (1,:) char = fullfile(wob.projectRoot, 'config', 'params_vidd.json')
end

if ~isfile(configPath)
    error('vidd:loadConfig:fileNotFound', 'Config file not found: %s', configPath);
end

raw = jsondecode(fileread(configPath));
cfg = raw;
cfg.configPath = configPath;
cfg.sources = struct();

for section = ["capacity", "disease"]
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

cfg.disease.d_disease_grid = reshape(raw.disease.d_disease_grid, 1, []);
end
