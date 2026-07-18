function cfg = loadConfig(configPath)
%LOADCONFIG Read the Model 1 configuration file into a struct.
%
%   cfg = wob.loadConfig() reads config/params_m1.json relative to the
%   project root.
%
%   cfg = wob.loadConfig(configPath) reads an explicit file.
%
%   The config is the single source of truth for every physiological and
%   device parameter. Scalar parameters under `patient` are stored as
%   {value, source} objects so that provenance travels with the number;
%   this function flattens them to plain scalars in cfg.patient while
%   preserving the provenance in cfg.sources.
%
%   See also wob.projectRoot

arguments
    configPath (1,:) char = fullfile(wob.projectRoot, 'config', 'params_m1.json')
end

if ~isfile(configPath)
    error('wob:loadConfig:fileNotFound', 'Config file not found: %s', configPath);
end

raw = jsondecode(fileread(configPath));

cfg = raw;
cfg.configPath = configPath;

% Flatten {value, source} wrappers in `patient` into plain numeric fields,
% keeping provenance in a parallel struct.
cfg.sources = struct();
patientFields = fieldnames(raw.patient);
for k = 1:numel(patientFields)
    name = patientFields{k};
    entry = raw.patient.(name);
    if isstruct(entry) && isfield(entry, 'value')
        cfg.patient.(name) = entry.value;
        cfg.sources.(name) = entry.source;
    end
end

% Device provenance travels the same way, resolved for the active
% coefficient set so that cfg.sources reflects what was actually used.
activeSet = raw.options_coefficient_set.active;
deviceNames = fieldnames(raw.devices);
for k = 1:numel(deviceNames)
    name = deviceNames{k};
    if isfield(raw.devices.(name), 'coefficients') && ...
       isfield(raw.devices.(name).coefficients, activeSet)
        cfg.sources.(name) = raw.devices.(name).coefficients.(activeSet).source;
    end
end

% jsondecode turns [1,2] into a column vector; the I:E ratio reads more
% naturally as a row.
cfg.pattern.IE_ratio = reshape(cfg.pattern.IE_ratio, 1, []);
end
