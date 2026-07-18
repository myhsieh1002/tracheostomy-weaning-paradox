function out = steadyState(cfgCO2, cfgM1, deviceName, opts)
%STEADYSTATE Steady-state PaCO2 at a prescribed ventilation (H1).
%
%   out = co2.steadyState(cfgCO2, cfgM1, deviceName) evaluates the alveolar
%   gas equation at the configured V_T and RR.
%
%   out = co2.steadyState(..., V_T=..., RR=..., VCO2=...) overrides.
%
%       V_A     = (V_T - V_D_total) * RR          L/min
%       PaCO2   = K * VCO2 / V_A                  mmHg,  K = 0.863
%
%   H1 ASKS A COUNTERFACTUAL, AND SHOULD SAY SO
%   -------------------------------------------
%   "At the same minute ventilation, the tracheostomy gives a lower PaCO2"
%   is true of this equation and is NOT what happens to a patient. Chadda
%   2002 measured both states in the same subjects: dead space fell 74 mL and
%   **PaCO2 did not change** -- the patients raised V_T from 330 to 400 mL
%   instead, at unchanged rate. PaCO2 is the regulated variable; ventilation
%   is the response.
%
%   So H1 is the mechanism isolated with the controller switched off, and
%   H2 (co2.requiredVentilation) is what the patient actually does. Both are
%   worth reporting; only the second is a clinical prediction.
%
%   Fields of `out`: PaCO2, V_A, V_E, V_T, RR, VCO2, ds.
%
%   See also co2.requiredVentilation, co2.deadSpaceTotal

arguments
    cfgCO2     (1,1) struct
    cfgM1      (1,1) struct
    deviceName (1,:) char
    opts.V_T  (1,1) double = NaN
    opts.RR   (1,1) double = NaN
    opts.VCO2 (1,1) double = NaN
    opts.Phi  (1,1) double = NaN
end

V_T  = pick(opts.V_T,  cfgCO2.pattern.V_T_L);
RR   = pick(opts.RR,   cfgCO2.pattern.RR);
VCO2 = pick(opts.VCO2, cfgCO2.metabolism.VCO2_mL_min);
K    = cfgCO2.constants.alveolar_K;

ds = co2.deadSpaceTotal(cfgCO2, cfgM1, deviceName, V_T, opts.Phi);

if V_T <= ds.total
    error('co2:steadyState:tidalBelowDeadSpace', ...
        ['V_T (%.3f L) does not exceed dead space (%.3f L) for device "%s": alveolar ' ...
         'ventilation is zero and PaCO2 is unbounded.'], V_T, ds.total, deviceName);
end

V_A = (V_T - ds.total) * RR;

out.PaCO2 = K * VCO2 / V_A;
out.V_A   = V_A;
out.V_E   = V_T * RR;
out.V_T   = V_T;
out.RR    = RR;
out.VCO2  = VCO2;
out.ds    = ds;
out.device = deviceName;
end

function v = pick(given, fallback)
if isnan(given), v = fallback; else, v = given; end
end
