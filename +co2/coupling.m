function cp = coupling(cfgCO2, cfgM1, cfgM2, deviceName, opts)
%COUPLING Model 1b -> Model 1 -> Model 2, end to end.
%
%   cp = co2.coupling(cfgCO2, cfgM1, cfgM2, deviceName) runs the full chain:
%   CO2 kinetics set the ventilation required to hold the target PaCO2,
%   Model 1 computes the effort of delivering it, and Model 2 receives the
%   result as its load.
%
%       required V_E  <- co2.requiredVentilation   (this model)
%       P_mus, WOB    <- wob.simulateEffort        (Model 1)
%       L             <- lc.utilisation            (Model 2)
%
%   Name-value: VCO2, VD_VT, RR, PaCO2_target, and patient overrides for
%   Model 1 (C_rs, R_aw_native, ...) via PatientOverrides.
%
%   WHY THE VENTILATION IS COMPUTED HERE AND NOT IN MODEL 1
%   -------------------------------------------------------
%   Model 1 solves V_T from a configured target_VA using SERIES dead space
%   alone. This model derives V_A from CO2 production and adds the alveolar
%   dead space of V/Q mismatch. Letting Model 1 re-solve would use a
%   different dead space than the one this model just computed -- two answers
%   to one question. So the ventilation is passed through explicitly.
%
%   WHAT THIS CLOSES
%   ----------------
%   Model 1 currently carries target_VA = 6.3 L/min as a grade-B derivation
%   from the VCO2 proportionality, and Sobol ranks it the single largest
%   driver of total load (ST = 0.56). This chain replaces that derivation
%   with a computed value, so the model's biggest lever stops being an
%   assumption. cp.target_VA_implied reports what Model 1 should be
%   configured with if run standalone.
%
%   Fields of `cp`: vent, effort, L_total, f_device, target_VA_implied,
%   PaCO2, VCO2, VD_VT, convention.
%
%   See also co2.requiredVentilation, wob.simulateEffort, lc.utilisation

arguments
    cfgCO2     (1,1) struct
    cfgM1      (1,1) struct
    cfgM2      (1,1) struct
    deviceName (1,:) char
    opts.VCO2              (1,1) double = NaN
    opts.VD_VT             (1,1) double = NaN
    opts.RR                (1,1) double = NaN
    opts.PaCO2_target      (1,1) double = NaN
    opts.PatientOverrides  (1,1) struct = struct()
end

if ~isnan(opts.VD_VT)
    cfgCO2.deadspace.VD_VT = opts.VD_VT;
end

% --- 1b: what ventilation does this patient need? ---
vent = co2.requiredVentilation(cfgCO2, cfgM1, deviceName, ...
    RR=opts.RR, VCO2=opts.VCO2, PaCO2_target=opts.PaCO2_target);

% --- Model 1: what does delivering it cost? ---
ventForM1 = struct('V_T', vent.V_T, 'RR', vent.RR, 'V_D', vent.V_D, ...
                   'V_A', vent.V_A, 'V_E', vent.V_E);
effort = wob.simulateEffort(cfgM1, deviceName, opts.PatientOverrides, ...
                            Ventilation=ventForM1);

% --- Model 2: the load, in its own convention ---
[L_total, meta] = lc.utilisation(effort, cfgM2.load_from_model1.use);

cp.device            = deviceName;
cp.vent              = vent;
cp.effort            = effort;
cp.L_total           = L_total;
cp.f_device          = effort.f_device;
cp.convention        = meta.convention;
cp.u_crit_pairing    = meta.u_crit_pairing;
cp.PaCO2             = vent.PaCO2;
cp.VCO2              = vent.VCO2;
cp.VD_VT             = cfgCO2.deadspace.VD_VT;
cp.target_VA_implied = vent.V_A;
end
