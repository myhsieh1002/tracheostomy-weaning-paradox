function cp = coupling(cfg1, cfg2, deviceName, patientOverrides, opts)
%COUPLING End-to-end Model 1 -> Model 2 load transfer.
%
%   cp = lc.coupling(cfg1, cfg2, deviceName) runs Model 1's Mode B for the
%   given device and returns the load decomposition that Model 2 consumes.
%
%   cp = lc.coupling(cfg1, cfg2, deviceName, patientOverrides) sweeps the
%   patient axis (e.g. struct('C_rs', 0.03, 'R_aw_native', 10)).
%
%   This is the only bridge between the two models. Nothing here hard-codes
%   a load: L comes from wob.simulateEffort every time, which is what makes
%   the two-model chain a genuine coupling rather than two separate studies
%   sharing a figure caption.
%
%   The device's share of the load is apportioned by the work fraction:
%
%       L_device  = f_device * L_total
%       L_disease = L_total - L_device
%
%   A caveat worth stating in the manuscript: f_device is a fraction of
%   WORK, while L is a PRESSURE. Splitting a pressure by a work fraction is
%   exact only if the pressure components are in fixed proportion across the
%   breath, which they are not -- the resistive terms peak at peak flow and
%   the elastic term at end-inspiration. cp.L_device_direct gives the
%   alternative, computed by re-running Mode B with the device's resistance
%   removed; the difference between the two is reported as
%   cp.apportionmentGap so the approximation is measured rather than assumed.
%
%   Fields of `cp`: L_total, L_device, L_disease, L_device_direct,
%   apportionmentGap, f_device, effort, meta, convention.
%
%   See also wob.simulateEffort, lc.utilisation, lc.rescueOutcome

%   The direct-apportionment diagnostic costs a second full Mode B solve, so
%   it is opt-in via Apportionment=true. Grid sweeps that need only L_total
%   should leave it off.

arguments
    cfg1             (1,1) struct
    cfg2             (1,1) struct
    deviceName       (1,:) char
    patientOverrides (1,1) struct = struct()
    opts.Apportionment (1,1) logical = false
end

useWhich = cfg2.load_from_model1.use;

effort = wob.simulateEffort(cfg1, deviceName, patientOverrides);
[L_total, meta] = lc.utilisation(effort, useWhich);

% Work-fraction apportionment (the spec's route).
L_device  = effort.f_device * L_total;
L_disease = L_total - L_device;

cp.device           = deviceName;
cp.L_total          = L_total;
cp.L_device         = L_device;
cp.L_disease        = L_disease;
cp.f_device         = effort.f_device;
cp.effort           = effort;
cp.convention       = meta.convention;
cp.u_crit_pairing   = meta.u_crit_pairing;
cp.meta             = meta;

if opts.Apportionment
    % Direct route: what the load would be with a zero-resistance device of
    % the same dead space. The difference is the device's actual pressure
    % contribution, without assuming proportionality.
    cfgNoDev = cfg1;
    activeSet = cfg1.options_coefficient_set.active;
    cfgNoDev.devices.(deviceName).coefficients.(activeSet).K1 = 0;
    cfgNoDev.devices.(deviceName).coefficients.(activeSet).K2 = 0;
    effortNoDev = wob.simulateEffort(cfgNoDev, deviceName, patientOverrides);
    L_noDevice = lc.utilisation(effortNoDev, useWhich);

    cp.L_device_direct  = L_total - L_noDevice;
    cp.apportionmentGap = L_device - cp.L_device_direct;
else
    cp.L_device_direct  = NaN;
    cp.apportionmentGap = NaN;
end
end
