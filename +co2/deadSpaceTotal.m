function ds = deadSpaceTotal(cfgCO2, cfgM1, deviceName, V_T, phi)
%DEADSPACETOTAL Physiologic dead space: series (from Model 1) + alveolar.
%
%   ds = co2.deadSpaceTotal(cfgCO2, cfgM1, deviceName, V_T) returns the dead
%   space budget in litres.
%
%   ds = co2.deadSpaceTotal(..., phi) reuses a precomputed alveolar fraction
%   rather than deriving it (the derivation reads the config and Model 1, so
%   hoisting it out of a sweep is worth it).
%
%       V_D_series   = anatomic - bypassed + apparatus     [from Model 1]
%       V_D_alveolar = phi * (V_T - V_D_series)            [V/Q mismatch]
%       V_D_total    = V_D_series + V_D_alveolar
%
%   THE SERIES TERM IS NOT RECOMPUTED HERE
%   --------------------------------------
%   It comes from wob.deadSpace(). Model 1b's spec (section 8) requires dV_D
%   to be single-sourced with Model 1; calling Model 1's function is the only
%   way that stays true as either model changes. It also inherits Model 1's
%   correction: an ETT and a tracheostomy tube bypass the SAME airway, so the
%   bypass term cancels between them and only apparatus volume differs.
%
%   WHY THE ALVEOLAR TERM SCALES WITH (V_T - V_D_series) AND NOT WITH V_T
%   ---------------------------------------------------------------------
%   See co2.alveolarDeadSpaceFraction. In short: making it a fraction of V_T
%   would absorb the device's volume and leave ventilatory demand completely
%   independent of the tube, which is false and would make H3 vacuous.
%
%   Fields of `ds` (litres): series, alveolar, total, VD_VT_eff, phi.
%
%   See also wob.deadSpace, co2.alveolarDeadSpaceFraction

arguments
    cfgCO2     (1,1) struct
    cfgM1      (1,1) struct
    deviceName (1,:) char
    V_T        (1,1) double {mustBePositive}
    phi        (1,1) double = NaN
end

if isnan(phi)
    phi = co2.alveolarDeadSpaceFraction(cfgCO2, cfgM1);
end

dev = wob.getDevice(cfgM1, deviceName);
m1  = wob.deadSpace(cfgM1, dev);
series = m1.total;

if V_T <= series
    error('co2:deadSpaceTotal:tidalBelowSeries', ...
        ['V_T (%.3f L) does not exceed the series dead space (%.3f L) for device "%s": ' ...
         'no gas reaches the alveoli.'], V_T, series, deviceName);
end

alveolar = phi * (V_T - series);

ds.series    = series;
ds.alveolar  = alveolar;
ds.total     = series + alveolar;
ds.VD_VT_eff = ds.total / V_T;
ds.phi       = phi;
ds.m1        = m1;
ds.device    = deviceName;
end
