function dev = getDevice(cfg, name, coefficientSet)
%GETDEVICE Fetch one airway device, resolving the active coefficient set.
%
%   dev = wob.getDevice(cfg, 'TRACH_8_0') returns a struct with fields
%   name, type, ID_mm, length_cm, K1, K2, Vd_apparatus_mL, grade, source.
%
%   dev = wob.getDevice(cfg, name, 'flevari_ett') overrides the active set
%   for a sensitivity run.
%
%   WHICH COEFFICIENTS
%   ------------------
%   Each device carries K1/K2 under more than one bench source, because the
%   ETT-vs-tracheostomy contrast is the scientific question and it must not
%   be contaminated by between-rig differences. The active set is
%   cfg.options_coefficient_set.active; 'guttmann_within_rig' is primary
%   because Guttmann 1993 measured both device classes on one rig, so any
%   rig-specific bias cancels from the contrast. See the _device_notes block
%   in config/params_m1.json.
%
%   See also wob.rohrerDrop, wob.loadConfig

arguments
    cfg            (1,1) struct
    name           (1,:) char
    coefficientSet (1,:) char = ''
end

if ~isfield(cfg.devices, name)
    error('wob:getDevice:unknownDevice', ...
        'Device "%s" is not defined in %s. Available: %s', ...
        name, cfg.configPath, strjoin(fieldnames(cfg.devices)', ', '));
end

entry = cfg.devices.(name);

if isempty(coefficientSet)
    coefficientSet = cfg.options_coefficient_set.active;
end

if ~isfield(entry.coefficients, coefficientSet)
    error('wob:getDevice:unknownCoefficientSet', ...
        'Coefficient set "%s" is not defined for device "%s". Available: %s', ...
        coefficientSet, name, strjoin(fieldnames(entry.coefficients)', ', '));
end

coef = entry.coefficients.(coefficientSet);

dev                 = rmfield(entry, 'coefficients');
dev.name            = name;
dev.K1              = coef.K1;
dev.K2              = coef.K2;
dev.length_cm       = coef.length_cm;
dev.grade           = coef.grade;
dev.source          = coef.source;
dev.coefficientSet  = coefficientSet;
end
