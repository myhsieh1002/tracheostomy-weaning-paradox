function out = requiredVentilation(cfgCO2, cfgM1, deviceName, opts)
%REQUIREDVENTILATION Ventilation needed to hold the target PaCO2 (H2).
%
%   out = co2.requiredVentilation(cfgCO2, cfgM1, deviceName) inverts the
%   alveolar gas equation for the target PaCO2 and solves for the tidal
%   volume that delivers it through this device's dead space.
%
%       V_A_req = K * VCO2 / PaCO2_target
%       V_A     = (1 - phi) * (V_T - V_D_series) * RR
%   =>  V_T     = V_A_req / ((1 - phi) * RR) + V_D_series
%       V_E     = V_T * RR
%
%   Closed form -- no iteration. The implicit-looking dependence (V_D depends
%   on V_T, V_T depends on V_D) dissolves once the alveolar dead space is
%   written as a fraction of alveolar gas rather than of tidal volume; see
%   co2.alveolarDeadSpaceFraction.
%
%   THIS IS THE PHYSIOLOGICALLY REAL QUESTION
%   -----------------------------------------
%   PaCO2 is the regulated variable, so ventilation is what moves. Chadda
%   2002 measured exactly this: removing 74 mL of dead space left PaCO2 and
%   RR unchanged and raised V_T from 330 to 400 mL. H1 (co2.steadyState) asks
%   the counterfactual with the controller switched off; this asks what the
%   patient actually does.
%
%   THE DEVICE ENTERS THROUGH V_D_series ONLY
%   -----------------------------------------
%   So the ventilatory saving of a tracheostomy over an ETT is
%
%       dV_E = (V_D_series,ETT - V_D_series,TRACH) * RR
%
%   i.e. the apparatus difference times the rate -- inheriting Model 1's
%   finding that the upper-airway bypass cancels between the two tubes. The
%   saving against the NON-INTUBATED state is much larger, and that is the
%   comparison the 72 mL literature actually describes.
%
%   Name-value: RR, VCO2, PaCO2_target, Phi override the config.
%
%   Fields of `out`: V_T, RR, V_E, V_A, V_D, VCO2, PaCO2, ds.
%
%   See also co2.steadyState, co2.coupling, co2.alveolarDeadSpaceFraction

arguments
    cfgCO2     (1,1) struct
    cfgM1      (1,1) struct
    deviceName (1,:) char
    opts.RR           (1,1) double = NaN
    opts.VCO2         (1,1) double = NaN
    opts.PaCO2_target (1,1) double = NaN
    opts.Phi          (1,1) double = NaN
end

RR     = pick(opts.RR,           cfgCO2.pattern.RR);
VCO2   = pick(opts.VCO2,         cfgCO2.metabolism.VCO2_mL_min);
target = pick(opts.PaCO2_target, cfgCO2.targets.PaCO2_target_mmHg);
K      = cfgCO2.constants.alveolar_K;

phi = opts.Phi;
if isnan(phi)
    phi = co2.alveolarDeadSpaceFraction(cfgCO2, cfgM1);
end

V_A_req = K * VCO2 / target;                       % L/min

dev        = wob.getDevice(cfgM1, deviceName);
V_D_series = wob.deadSpace(cfgM1, dev).total;      % L

V_T = V_A_req / ((1 - phi) * RR) + V_D_series;

ds = co2.deadSpaceTotal(cfgCO2, cfgM1, deviceName, V_T, phi);

out.V_T    = V_T;
out.RR     = RR;
out.V_E    = V_T * RR;
out.V_A    = (V_T - ds.total) * RR;
out.V_D    = ds.total;
out.VCO2   = VCO2;
out.PaCO2  = target;
out.phi    = phi;
out.ds     = ds;
out.device = deviceName;
end

function v = pick(given, fallback)
if isnan(given), v = fallback; else, v = given; end
end
