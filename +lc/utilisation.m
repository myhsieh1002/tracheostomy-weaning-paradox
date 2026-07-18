function [L, meta] = utilisation(effort, useWhich)
%UTILISATION Extract the load L from a Model 1 effort result, safely paired.
%
%   [L, meta] = lc.utilisation(effort, useWhich) pulls the requested load
%   convention out of a wob.simulateEffort result and returns, alongside it,
%   the fatigue threshold that MUST be used with it.
%
%   The pairing is the whole point of this function. u = L/C is
%   dimensionless, so nothing in the arithmetic stops a duty-weighted load
%   being compared against a threshold measured for a mean inspiratory
%   pressure. The two conventions differ by the duty cycle -- roughly a
%   factor of three -- so mixing them silently mis-scales the entire
%   bifurcation.
%
%       useWhich          L is...                    pair with u_crit =
%       ---------------   -------------------------  ------------------
%       'mean'            mean INSPIRATORY P_mus     0.40  Roussos & Macklem 1977
%                         => u = Pi/Pimax                  (PMID 893274)
%       'duty_weighted'   PTP/Ttot                   0.15  Bellemare & Grassino 1982
%                         => u = tension-time index        (PMID 7174413)
%       'peak'            peak P_mus                 NO PUBLISHED THRESHOLD
%
%   'peak' is accepted for exploration but carries no anchored threshold and
%   returns NaN for u_crit; callers must supply one explicitly and say so.
%
%   Fields of `meta`: convention, u_crit_pairing, source, note.
%
%   See also wob.simulateEffort, lc.coupling

arguments
    effort   (1,1) struct
    useWhich (1,:) char
end

switch lower(useWhich)
    case 'mean'
        L = effort.P_mus_mean;
        meta.convention     = "mean inspiratory P_mus";
        meta.u_crit_pairing = 0.40;
        meta.source         = "Roussos & Macklem 1977, PMID 893274 (Pdicrit, a pure pressure ratio)";
        meta.note           = "u = L/C_max is Pi/Pimax.";
    case {'duty_weighted', 'dutyweighted'}
        L = effort.P_mus_dutyWeighted;
        meta.convention     = "duty-weighted P_mus (PTP/Ttot)";
        meta.u_crit_pairing = 0.15;
        meta.source         = "Bellemare & Grassino 1982, PMID 7174413 (TTdi)";
        meta.note           = "u = L/C_max is the tension-time index, which already carries Ti/Ttot.";
    case 'peak'
        L = effort.P_mus_peak;
        meta.convention     = "peak inspiratory P_mus";
        meta.u_crit_pairing = NaN;
        meta.source         = "none - no published threshold exists for a peak-pressure ratio";
        meta.note           = "Exploration only. Supply u_crit explicitly and justify it.";
    otherwise
        error('lc:utilisation:unknownConvention', ...
            'load_from_model1.use must be ''mean'', ''duty_weighted'' or ''peak'' (got "%s").', useWhich);
end
end
