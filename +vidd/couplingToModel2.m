function cp = couplingToModel2(cfgV, cfgM1, cfgM2, scenarioName, opts)
%COUPLINGTOMODEL2 The dynamic rescue window: C_max(t) -> Model 2.
%
%   cp = vidd.couplingToModel2(cfgV, cfgM1, cfgM2, scenarioName) simulates
%   the capacity trajectory under a ventilation scenario and walks it
%   through Model 2's fold, producing the rescue window as a function of
%   time (Model 2's m2_exp6 / this model's vidd_exp4 -- the same figure).
%
%   THE DEVICE MOVES BOTH AXES, AND THE SPEC ONLY SAYS SO IN PASSING
%   ---------------------------------------------------------------
%   A tracheostomy acts on this programme's two axes at once:
%
%     LOAD     via Model 1 -- a shorter, wider tube costs less pressure.
%     CAPACITY via this model -- it enables lighter sedation, which moves
%              diaphragm activity A towards its optimum and slows the
%              degradation term.
%
%   So the scenario carries BOTH a device (which sets L through Model 1) and
%   a support level (which sets A here). Both boundaries of the rescue window
%   move: the fold falls as C_max(t) falls, and the operating point moves as
%   L changes. Neither can be held fixed while asking about the other, which
%   is why this function takes a scenario rather than a device.
%
%   THE SLOW-FAST SEPARATION IS WHAT MAKES THIS LEGITIMATE
%   -----------------------------------------------------
%   Capacity moves on DAYS (tau ~ 1-2 d here); Model 2's fatigue dynamics
%   move on HOURS (tau = 1/alpha ~ 7 h). The ratio is ~5-50x, so at each
%   C_max(t) Model 2's fast dynamics have effectively equilibrated and its
%   fold can be recomputed quasi-statically. cp.separation reports the
%   actual ratio so the assumption is measured, not asserted; below ~5 the
%   quasi-static reading stops being safe and the two would need solving
%   together.
%
%   THE FOLD IS NOT RE-DERIVED HERE
%   -------------------------------
%   It comes from lc.rescueWindowTrajectory, which calls lc.rescueWindow --
%   the named API in Model 2's spec section 2.5b. Model 2's section 9 is
%   explicit that 2b must not write its own fold solve.
%
%   Name-value: Days, C0, PatientOverrides (Model 1 patient axes).
%
%   Fields of `cp`: traj (the Model 2 trajectory), cap (the capacity
%   trajectory), L_ETT, L_TRACH, scenario, separation.
%
%   See also lc.rescueWindowTrajectory, vidd.simulateCapacity, lc.coupling

arguments
    cfgV          (1,1) struct
    cfgM1         (1,1) struct
    cfgM2         (1,1) struct
    scenarioName  (1,:) char
    opts.Days             (1,1) double = NaN
    opts.C0               (1,1) double = NaN
    opts.PatientOverrides (1,1) struct = struct()
    opts.Params           (1,1) struct = struct()
end

if ~isfield(cfgV.strategy.scenarios, scenarioName)
    error('vidd:couplingToModel2:unknownScenario', ...
        'Unknown scenario "%s". Available: %s', scenarioName, ...
        strjoin(fieldnames(cfgV.strategy.scenarios)', ', '));
end

scenario = cfgV.strategy.scenarios.(scenarioName);

% Calibrate, THEN apply overrides -- vidd.params alone leaves k_syn/g_min/
% p_mono as NaN (they are solved by vidd.calibrate), which would make the
% whole capacity trajectory NaN and fail downstream positivity checks.
p = vidd.calibrate(cfgV);
fn = fieldnames(opts.Params);
for k = 1:numel(fn)
    p.(fn{k}) = opts.Params.(fn{k});
end

pM2 = lc.calibrate(cfgM2);

% --- capacity axis: this model ---
A = vidd.supportToActivity(scenario.support_level, cfgV);
cap = vidd.simulateCapacity(cfgV, p, @(t) A, C0=opts.C0, Days=opts.Days);

% --- load axis: Model 1, never hard-coded ---
devETT   = cfgV.coupling.devices{1};
devTRACH = cfgV.coupling.devices{2};
L_ETT   = lc.coupling(cfgM1, cfgM2, devETT,   opts.PatientOverrides).L_total;
L_TRACH = lc.coupling(cfgM1, cfgM2, devTRACH, opts.PatientOverrides).L_total;

% --- the seam: Model 2's named API, once per time point ---
traj = lc.rescueWindowTrajectory(cap.t, cap.C, [L_TRACH, L_ETT], pM2);

% --- is the quasi-static reading actually safe here? ---
tau_slow_days = mean(cap.tau_inst);
tau_fast_days = (1 / pM2.alpha) / 24;
cp.separation = tau_slow_days / tau_fast_days;
if cp.separation < 5
    warning('vidd:couplingToModel2:weakSeparation', ...
        ['Timescale separation is only %.1fx (capacity tau %.2f d vs Model 2 tau %.2f d). ' ...
         'The quasi-static coupling assumes the fast system equilibrates between capacity ' ...
         'steps; below ~5x that reading is not safe and the two should be solved together.'], ...
        cp.separation, tau_slow_days, tau_fast_days);
end

cp.traj          = traj;
cp.cap           = cap;
cp.L_ETT         = L_ETT;
cp.L_TRACH       = L_TRACH;
cp.dL_device     = L_ETT - L_TRACH;
cp.A             = A;
cp.scenario      = scenario;
cp.scenarioName  = scenarioName;
cp.p             = p;
cp.pM2           = pM2;
cp.tau_slow_days = tau_slow_days;
cp.tau_fast_days = tau_fast_days;
end
