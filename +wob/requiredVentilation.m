function vent = requiredVentilation(cfg, dev, pattern)
%REQUIREDVENTILATION Ventilation needed to hit the target alveolar ventilation.
%
%   vent = wob.requiredVentilation(cfg, dev, pattern) inverts
%
%       V_A = (V_T - V_D) * RR
%
%   for the breathing pattern that delivers cfg.patient.target_VA_L_min
%   with the dead space imposed by `dev`.
%
%   Which variable is solved for is set by cfg.options.mode_B_fix:
%       'RR'   hold RR at pattern.RR and solve for V_T   (default)
%       'V_T'  hold V_T at pattern.V_T_L and solve for RR
%
%   Holding RR and solving V_T is the default because respiratory rate is
%   the more tightly regulated of the two in a weaning patient, and because
%   it keeps the inspiratory time fixed across devices so that the device
%   contrast is not confounded by a changed duty cycle.
%
%   Fields of `vent`: V_T, RR, V_E, V_A, V_D, pattern (updated).
%
%   See also wob.deadSpace, wob.simulateModeB

arguments
    cfg     (1,1) struct
    dev     (1,1) struct
    pattern (1,1) struct = cfg.pattern
end

ds = wob.deadSpace(cfg, dev);
targetVA = cfg.patient.target_VA_L_min;

switch upper(string(cfg.options.mode_B_fix))
    case "RR"
        RR  = pattern.RR;
        V_T = targetVA / RR + ds.total;
    case {"V_T", "VT"}
        V_T = pattern.V_T_L;
        if V_T <= ds.total
            error('wob:requiredVentilation:tidalBelowDeadSpace', ...
                ['V_T (%.3f L) does not exceed dead space (%.3f L) for device "%s": ' ...
                 'no alveolar ventilation is possible at any rate.'], V_T, ds.total, dev.name);
        end
        RR = targetVA / (V_T - ds.total);
    otherwise
        error('wob:requiredVentilation:badFixOption', ...
            'cfg.options.mode_B_fix must be ''RR'' or ''V_T'' (got "%s").', ...
            cfg.options.mode_B_fix);
end

pattern.V_T_L = V_T;
pattern.RR    = RR;

vent.V_T     = V_T;
vent.RR      = RR;
vent.V_E     = V_T * RR;
vent.V_A     = (V_T - ds.total) * RR;
vent.V_D     = ds.total;
vent.ds      = ds;
vent.pattern = pattern;
end
