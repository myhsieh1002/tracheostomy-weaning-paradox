function res = simulateModeB(cfg, dev, pattern, opts)
%SIMULATEMODEB Inverse-solve the muscle pressure needed to sustain target V_A.
%
%   res = wob.simulateModeB(cfg, dev, pattern) is the active-drive mode
%   behind H3 and H4, and the sole supplier of load to Model 2.
%
%   METHOD
%   ------
%   The device's dead space fixes the ventilation required to reach the
%   target alveolar ventilation (wob.requiredVentilation). That pattern
%   prescribes flow and volume, and the equation of motion is then inverted
%   ALGEBRAICALLY for the driving pressure:
%
%       P_mus(t) + P_vent(t) = R_total(V')*V' + V/C + I*V'' + P0
%
%   With P_vent = 0 (T-piece / trach mask, the weaning condition), every
%   term on the right is known once the waveform is prescribed, so
%
%       P_mus(t) = R_aw*V' + dP_device(V') + V/C + I*V'' + P0
%
%   No ODE integration is needed. The build spec anticipated an ODE solve
%   here, but prescribing the waveform makes the equation of motion
%   explicit in P_mus: it is a direct evaluation, exact to machine
%   precision, with no solver tolerance to tune. The ODE form would be
%   required only for the forward problem (given P_mus, find V), which this
%   model never poses.
%
%   THE DECOMPOSITION
%   -----------------
%   Each pressure term is integrated against the same flow with the same
%   quadrature rule, so the work components sum to the total by
%   construction:
%
%       WOB_total = WOB_elastic + WOB_resistive_native + WOB_device
%                   (+ WOB_inertial, when enabled)
%
%   f_device = WOB_device / WOB_total is the H4 metric.
%
%   Fields of `res` include the pattern, the pressure traces, the work
%   decomposition in J, J/L and J/min, the PTP, the P_mus summaries
%   (peak / mean / duty-weighted mean) and f_device.
%
%   EXTERNALLY SUPPLIED VENTILATION
%   -------------------------------
%   res = wob.simulateModeB(cfg, dev, pattern, Ventilation=v) skips the
%   internal wob.requiredVentilation and uses the supplied struct (fields
%   V_T, RR, V_D, V_A, V_E) instead.
%
%   This exists for Model 1b. Model 1 knows only SERIES dead space (anatomic
%   + apparatus); Model 1b adds the alveolar dead space of V/Q mismatch and
%   derives the required ventilation from CO2 kinetics rather than from a
%   configured target_VA. If Model 1 re-solved the ventilation itself, the
%   two models would silently disagree about the dead space -- which is
%   exactly what the Model 1b spec's section 8 forbids. So the division of
%   labour is: 1b computes the ventilation, Model 1 computes the mechanics
%   of delivering it.
%
%   See also wob.simulateModeA, wob.requiredVentilation, wob.simulateEffort,
%   co2.coupling

arguments
    cfg     (1,1) struct
    dev     (1,1) struct
    pattern (1,1) struct = cfg.pattern
    opts.Ventilation (1,1) struct = struct()
end

c = wob.constants();

if isempty(fieldnames(opts.Ventilation))
    vent = wob.requiredVentilation(cfg, dev, pattern);
else
    vent = opts.Ventilation;
    required = ["V_T","RR","V_D","V_A","V_E"];
    missing = required(~isfield(vent, required));
    if ~isempty(missing)
        error('wob:simulateModeB:incompleteVentilation', ...
            'Supplied Ventilation is missing field(s): %s', strjoin(missing, ', '));
    end
    vent.pattern       = pattern;
    vent.pattern.V_T_L = vent.V_T;
    vent.pattern.RR    = vent.RR;
end

pat = wob.breathingPattern(vent.pattern);

C_rs  = cfg.patient.C_rs;
R_aw  = cfg.patient.R_aw_native;
I     = cfg.patient.inertance;
P0    = cfg.patient.P0;

if C_rs <= 0
    error('wob:simulateModeB:badCompliance', 'C_rs must be positive (got %g).', C_rs);
end

% --- Pressure components over the whole cycle ---
P_elastic         = pat.volume / C_rs;
P_resistiveNative = R_aw * pat.flow;
P_device          = wob.rohrerDrop(dev, pat.flow);

if cfg.options.use_inertance
    P_inertial = I * pat.accel;
else
    P_inertial = zeros(size(pat.flow));
end

P_mus = P_resistiveNative + P_device + P_elastic + P_inertial + P0;

insp = pat.isInsp;

% --- Work decomposition (cmH2O*L -> J) ---
toJ = c.CMH2O_L_TO_JOULE;

W_elastic  = wob.workIntegral(pat.t, P_elastic,         pat.flow, insp) * toJ;
W_native   = wob.workIntegral(pat.t, P_resistiveNative, pat.flow, insp) * toJ;
W_device   = wob.workIntegral(pat.t, P_device,          pat.flow, insp) * toJ;
W_inertial = wob.workIntegral(pat.t, P_inertial,        pat.flow, insp) * toJ;
W_offset   = wob.workIntegral(pat.t, P0*ones(size(pat.flow)), pat.flow, insp) * toJ;
W_total    = wob.workIntegral(pat.t, P_mus,             pat.flow, insp) * toJ;

% --- PTP over inspiration ---
PTP_breath = wob.ptpIntegral(pat.t, P_mus, insp);

% --- P_mus summaries ---
% `mean` is the mean over INSPIRATION and is the quantity that pairs with
% MIP to give u = L/C = Pi/Pimax, the ratio for which u_crit ~ 0.4 was
% measured (Roussos & Macklem 1977). `dutyWeighted` divides the same
% integral by the whole cycle, reproducing the tension-time index for which
% the threshold is ~0.15 (Bellemare & Grassino 1982). The two are
% consistent parameterisations of the same physiology, related by the duty
% cycle; Model 2 must pair each with its own threshold.
P_mus_peak         = max(P_mus(insp));
P_mus_mean         = PTP_breath / pat.Ti;
P_mus_dutyWeighted = PTP_breath / pat.Ttot;

res.device  = dev.name;
res.pat     = pat;
res.vent    = vent;

res.P_mus             = P_mus;
res.P_elastic         = P_elastic;
res.P_resistiveNative = P_resistiveNative;
res.P_device          = P_device;
res.P_inertial        = P_inertial;

res.WOB_total_J    = W_total;
res.WOB_elastic_J  = W_elastic;
res.WOB_native_J   = W_native;
res.WOB_device_J   = W_device;
res.WOB_inertial_J = W_inertial;
res.WOB_offset_J   = W_offset;

res.WOB_total_J_L    = W_total / pat.V_T;
res.WOB_total_J_min  = W_total * pat.RR;
res.WOB_device_J_L   = W_device / pat.V_T;
res.WOB_device_J_min = W_device * pat.RR;

res.f_device = W_device / W_total;

res.PTP_per_breath = PTP_breath;
res.PTP_per_min    = PTP_breath * pat.RR;

res.P_mus_peak         = P_mus_peak;
res.P_mus_mean         = P_mus_mean;
res.P_mus_dutyWeighted = P_mus_dutyWeighted;

res.V_T = pat.V_T;
res.RR  = pat.RR;
res.V_E = vent.V_E;
res.V_A = vent.V_A;
res.V_D = vent.V_D;
end
