function res = simulateModeA(cfg, dev, pattern)
%SIMULATEMODEA Imposed device load under a prescribed flow waveform.
%
%   res = wob.simulateModeA(cfg, dev, pattern) prescribes a physiological
%   flow waveform and integrates ONLY the pressure drop across the airway
%   device. This is the clean passive comparison behind H1 and H2: with
%   V_T and RR matched across devices, any difference in the result is
%   attributable to the tube alone.
%
%   Mode A deliberately says nothing about the patient's compliance or
%   native airway resistance -- it is not asked to. Use wob.simulateModeB
%   when the question involves total load or the device's share of it.
%
%   pattern defaults to cfg.pattern if omitted.
%
%   Fields of `res`:
%       pat            the breathing pattern used
%       dP_device      pressure drop across the device over the cycle (cmH2O)
%       WOB_device_J   imposed work per breath (J)
%       WOB_device_J_L imposed work per litre ventilated (J/L)
%       WOB_device_J_min imposed work per minute (J/min)
%       PTP_device_per_breath, PTP_device_per_min  (cmH2O*s, cmH2O*s/min)
%       peak_dP_device, mean_dP_device  over inspiration (cmH2O)
%       peak_flow      (L/s)
%
%   See also wob.simulateModeB, wob.rohrerDrop

arguments
    cfg     (1,1) struct
    dev     (1,1) struct
    pattern (1,1) struct = cfg.pattern
end

c = wob.constants();
pat = wob.breathingPattern(pattern);

dP_device = wob.rohrerDrop(dev, pat.flow);

insp = pat.isInsp;

% Work is integrated over inspiration only: on a T-piece, expiration is
% passive and its device pressure drop is paid for by stored elastic
% recoil, not by the inspiratory muscles.
W_cmH2O_L = wob.workIntegral(pat.t, dP_device, pat.flow, insp);
PTP_breath = wob.ptpIntegral(pat.t, dP_device, insp);

res.pat              = pat;
res.device           = dev.name;
res.dP_device        = dP_device;

res.WOB_device_J     = W_cmH2O_L * c.CMH2O_L_TO_JOULE;
res.WOB_device_J_L   = res.WOB_device_J / pat.V_T;
res.WOB_device_J_min = res.WOB_device_J * pat.RR;

res.PTP_device_per_breath = PTP_breath;
res.PTP_device_per_min    = PTP_breath * pat.RR;

res.peak_dP_device = max(dP_device(insp));
res.mean_dP_device = PTP_breath / pat.Ti;
res.peak_flow      = max(pat.flow(insp));

res.V_T = pat.V_T;
res.RR  = pat.RR;
res.V_E = pat.V_T * pat.RR;
end
