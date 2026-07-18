function phi = alveolarDeadSpaceFraction(cfgCO2, cfgM1)
%ALVEOLARDEADSPACEFRACTION Convert a conventional V_D/V_T into the model's phi.
%
%   phi = co2.alveolarDeadSpaceFraction(cfgCO2, cfgM1) returns the fraction
%   of ALVEOLAR gas that is ventilated but not perfused.
%
%   WHY phi AND NOT V_D/V_T DIRECTLY
%   --------------------------------
%   The conventional Bohr-Enghoff V_D/V_T is a fraction of TIDAL volume, and
%   it lumps two physically different things: the series dead space of tube
%   and conducting airways (a fixed volume, set by the device) and the
%   alveolar dead space of V/Q mismatch (set by the lung).
%
%   Modelling the alveolar part as a fraction of V_T destroys the result.
%   With V_D_total = VD_VT * V_T, the required ventilation is
%
%       V_E = V_A_req / (1 - VD_VT)
%
%   which contains NO device term at all: scaling the dead space with V_T
%   absorbs the tube's volume entirely, so swapping an ETT for a
%   tracheostomy would have exactly zero effect on ventilatory demand. That
%   is obviously false, and it would make H3 vacuous.
%
%   The physical parametrisation makes the alveolar part a fraction of the
%   ALVEOLAR gas -- the tidal volume that actually reaches alveoli:
%
%       V_D_alveolar = phi * (V_T - V_D_series)
%       V_A          = (1 - phi) * (V_T - V_D_series) * RR
%
%   Now phi is a property of the lung's V/Q distribution and V_D_series is a
%   property of the device, and the two are orthogonal -- which is what lets
%   the model ask a device question and a disease question separately.
%
%   THE CONVERSION
%   --------------
%   phi is calibrated so that the configured V_D/V_T is reproduced at a
%   reference condition (the ETT arm at the configured V_T):
%
%       VD_VT = [V_D_series + phi*(V_T - V_D_series)] / V_T
%   =>  phi   = (VD_VT*V_T - V_D_series) / (V_T - V_D_series)
%
%   A configured V_D/V_T below what the series dead space alone implies is
%   physically impossible -- gas must traverse the tube whatever V/Q does --
%   and is rejected rather than clamped, because silently clamping would
%   report a V_D/V_T the config did not ask for.
%
%   See also co2.deadSpaceTotal, co2.requiredVentilation

arguments
    cfgCO2 (1,1) struct
    cfgM1  (1,1) struct
end

V_T_ref = cfgCO2.pattern.V_T_L;
VD_VT   = cfgCO2.deadspace.VD_VT;
refDev  = cfgCO2.arms.ett;

dev    = wob.getDevice(cfgM1, refDev);
series = wob.deadSpace(cfgM1, dev).total;

if V_T_ref <= series
    error('co2:alveolarDeadSpaceFraction:tidalBelowSeries', ...
        'Reference V_T (%.3f L) does not exceed the series dead space (%.3f L).', ...
        V_T_ref, series);
end

phi = (VD_VT * V_T_ref - series) / (V_T_ref - series);

if phi < 0
    error('co2:alveolarDeadSpaceFraction:impossibleVDVT', ...
        ['Configured V_D/V_T = %.3f is below the %.3f implied by the series dead space ' ...
         'alone (%.1f mL at V_T = %.3f L, device %s). No V/Q distribution can achieve it: ' ...
         'the gas still has to traverse the tube and conducting airways.'], ...
        VD_VT, series/V_T_ref, series*1e3, V_T_ref, refDev);
end
if phi >= 1
    error('co2:alveolarDeadSpaceFraction:noAlveolarVentilation', ...
        ['Configured V_D/V_T = %.3f implies an alveolar dead-space fraction of %.3f: ' ...
         'no alveolar ventilation is possible at any tidal volume.'], VD_VT, phi);
end
end
