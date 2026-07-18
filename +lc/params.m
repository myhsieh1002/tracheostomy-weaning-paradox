function p = params(cfg, overrides)
%PARAMS Assemble the dynamical parameter struct, with optional overrides.
%
%   p = lc.params(cfg) returns struct with fields alpha, beta, s, u_crit.
%
%   p = lc.params(cfg, struct('beta', 1.5)) overrides selected parameters
%   without mutating cfg. Used by the grids and the sensitivity analysis.

arguments
    cfg       (1,1) struct
    overrides (1,1) struct = struct()
end

p = struct('alpha',  cfg.dynamics.alpha, ...
           'beta',   cfg.dynamics.beta, ...
           's',      cfg.dynamics.s, ...
           'u_crit', cfg.dynamics.u_crit);

names = fieldnames(overrides);
for k = 1:numel(names)
    name = names{k};
    if ~isfield(p, name)
        error('lc:params:unknownParameter', ...
            '"%s" is not a dynamical parameter. Valid: alpha, beta, s, u_crit.', name);
    end
    p.(name) = overrides.(name);
end

mustBePositive(p.alpha);
mustBePositive(p.beta);
mustBePositive(p.s);
mustBePositive(p.u_crit);
end
