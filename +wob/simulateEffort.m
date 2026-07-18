function out = simulateEffort(cfg, deviceName, overrides, opts)
%SIMULATEEFFORT Stable Model 1 -> Model 2 interface.
%
%   out = wob.simulateEffort(cfg, deviceName) runs Mode B and returns the
%   compact summary that Model 2 consumes as its load parameter.
%
%   out = wob.simulateEffort(cfg, deviceName, overrides) applies a struct
%   of patient-parameter overrides (e.g. struct('C_rs', 0.03,
%   'R_aw_native', 10)) before simulating, without mutating cfg. This is
%   how the disease-severity grid and the sensitivity analysis sweep the
%   patient axis.
%
%   Fields of `out`:
%       P_mus_peak, P_mus_mean, P_mus_dutyWeighted   (cmH2O)
%       f_device                                     (dimensionless)
%       WOB_total, WOB_device                        (J/breath)
%       WOB_total_J_min, WOB_device_J_min            (J/min)
%       PTP_per_min                                  (cmH2O*s/min)
%       V_E, V_T, RR, V_D
%
%   This function is the contract with Model 2: lc.coupling calls it and
%   nothing else in Model 1. Keeping the surface this narrow is what lets
%   either model be reworked without touching the other.
%
%   See also wob.simulateModeB, lc.coupling

%   out = wob.simulateEffort(..., Ventilation=v) passes an externally
%   computed ventilation through to Mode B instead of solving for it. This is
%   Model 1b's entry point -- see wob.simulateModeB for why the ventilation
%   must not be solved twice.

arguments
    cfg        (1,1) struct
    deviceName (1,:) char
    overrides  (1,1) struct = struct()
    opts.Ventilation (1,1) struct = struct()
end

% Apply overrides to a local copy so callers can sweep parameters freely.
fields = fieldnames(overrides);
for k = 1:numel(fields)
    name = fields{k};
    if ~isfield(cfg.patient, name)
        error('wob:simulateEffort:unknownOverride', ...
            ['"%s" is not a patient parameter. Valid: %s'], ...
            name, strjoin(fieldnames(cfg.patient)', ', '));
    end
    cfg.patient.(name) = overrides.(name);
end

dev = wob.getDevice(cfg, deviceName);
res = wob.simulateModeB(cfg, dev, cfg.pattern, Ventilation=opts.Ventilation);

out.device             = deviceName;
out.P_mus_peak         = res.P_mus_peak;
out.P_mus_mean         = res.P_mus_mean;
out.P_mus_dutyWeighted = res.P_mus_dutyWeighted;
out.f_device           = res.f_device;
out.WOB_total          = res.WOB_total_J;
out.WOB_device         = res.WOB_device_J;
out.WOB_elastic        = res.WOB_elastic_J;
out.WOB_native         = res.WOB_native_J;
out.WOB_total_J_min    = res.WOB_total_J_min;
out.WOB_device_J_min   = res.WOB_device_J_min;
out.PTP_per_min        = res.PTP_per_min;
out.V_E                = res.V_E;
out.V_T                = res.V_T;
out.RR                 = res.RR;
out.V_D                = res.V_D;
out.dutyCycle          = res.pat.dutyCycle;
end
