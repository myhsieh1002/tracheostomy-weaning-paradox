function p = params(cfg, overrides)
%PARAMS Assemble the VIDD parameter struct, with optional overrides.
%
%   p = vidd.params(cfg) returns the raw parameter struct. The values that
%   vidd.calibrate solves for (k_syn, g_min, p_mono) come back as NaN until
%   calibration runs -- use vidd.calibrate for a usable parameter set.
%
%   p = vidd.params(cfg, struct('d_disease', 0.1)) overrides without
%   mutating cfg.
%
%   Note p_exp/q_exp rather than p/q: `p` is the struct itself, and p.p
%   would be legal but unreadable.

arguments
    cfg       (1,1) struct
    overrides (1,1) struct = struct()
end

p = struct( ...
    'C_max0',    cfg.capacity.C_max0, ...
    'k_deg',     cfg.capacity.k_deg, ...
    'k_syn',     nanIfEmpty(cfg.capacity.k_syn), ...
    'h0',        cfg.capacity.h0, ...
    'g0',        cfg.capacity.g0, ...
    'g_mode',    string(cfg.modes.g_mode), ...
    'g_min',     nanIfEmpty(cfg.capacity.g_min), ...
    'p_mono',    nanIfEmpty(cfg.capacity.p_mono), ...
    'A_star',    cfg.capacity.A_star, ...
    'a_disuse',  cfg.capacity.a_disuse, ...
    'a_injury',  cfg.capacity.a_injury, ...
    'p_exp',     cfg.capacity.p, ...
    'q_exp',     cfg.capacity.q, ...
    'd_mode',    string(cfg.modes.d_mode), ...
    'd_disease', cfg.disease.d_disease, ...
    'sepsis_offset_fraction', cfg.disease.sepsis_offset_fraction);

names = fieldnames(overrides);
for k = 1:numel(names)
    name = names{k};
    if ~isfield(p, name)
        error('vidd:params:unknownParameter', ...
            '"%s" is not a VIDD parameter. Valid: %s', name, strjoin(fieldnames(p)', ', '));
    end
    p.(name) = overrides.(name);
end

mustBePositive(p.k_deg);
mustBeNonnegative(p.d_disease);
mustBeInRange(p.h0, 0, 1);
mustBeInRange(p.A_star, 0, 1);
mustBeInRange(p.sepsis_offset_fraction, 0, 1);

if ~ismember(p.g_mode, ["monotonic","ushape"])
    error('vidd:params:unknownGMode', ...
        'g_mode must be ''monotonic'' or ''ushape'' (got "%s").', p.g_mode);
end
if ~ismember(p.d_mode, ["offset","rate"])
    error('vidd:params:unknownDMode', ...
        'd_mode must be ''offset'' or ''rate'' (got "%s").', p.d_mode);
end
end

function v = nanIfEmpty(x)
% jsondecode turns JSON null into []; these are the fields vidd.calibrate
% fills in, and NaN says "not yet solved" more loudly than [] does.
if isempty(x), v = NaN; else, v = x; end
end
