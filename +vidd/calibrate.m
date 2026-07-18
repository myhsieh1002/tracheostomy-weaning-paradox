function [p, info] = calibrate(cfg, opts)
%CALIBRATE Solve k_syn, g_min and p_mono against published physiology.
%
%   [p, info] = vidd.calibrate(cfg) returns a usable parameter set.
%
%   WHAT IS ANCHORED, AND WHAT IS SOLVED
%   ------------------------------------
%   ANCHORED (not touched here):
%     k_deg = 0.064 /day -- Jaber 2010 (PMID 20813887): twitch airway
%       pressure fell 32+-6% over 6 days of ventilation, k = -ln(0.68)/6.
%       Force, in cmH2O, over days: the model's own state and timescale.
%       Corroborated by four independent thickness studies (0.048-0.115).
%       g(0) = 1 by convention, so k_deg IS the complete-disuse rate.
%
%   SOLVED HERE:
%     g_min, p_mono -- least squares against Zambon 2016's measured atrophy
%       rates by ventilatory mode, the only data stratifying RATE by support.
%     k_syn -- from Lecronier 2022's TWO groups (see below). Data, not a
%       modelling requirement.
%
%   k_syn IS IDENTIFIABLE AFTER ALL -- FROM TWO GROUPS, NOT ONE
%   ----------------------------------------------------------
%   Every positive rate in this literature is a NET rate, and a net rate
%   cannot separate k_syn from k_deg*g. That is true of any SINGLE group,
%   and it is why Zambon's +2.3%/day at spontaneous breathing is unusable
%   here.
%
%   Lecronier 2022 reports TWO groups re-measured at the same intervention
%   -- the switch to pressure support -- from different starting points:
%
%       septic:     Ptr,stim 6.3 cmH2O, +19% over 4 days
%       non-septic: Ptr,stim 9.8 cmH2O,  -7% over 4 days
%
%   Same activity means the same C* and the same tau, so two exponentials
%   through two points give two equations in two unknowns:
%
%       1.19*6.3 = C* + (6.3 - C*)*E
%       0.93*9.8 = C* + (9.8 - C*)*E,       E = exp(-4/tau)
%
%   Subtracting eliminates C*: E = 0.462, tau = 5.18 d, C* = 8.52 cmH2O.
%   With tau = 1/(syn + loss) and loss = k_deg*g(A_PS) already known, syn
%   follows, and k_syn = syn/h(A_PS).
%
%   The solve pays for itself twice: it also implies a septic starting
%   capacity of 65.5% of baseline -- a 34.5% offset -- against the 36% that
%   Demoule's independent regression coefficient gives. Two studies, two
%   methods, agreeing to two percentage points.
%
%   WHY LEVINE IS NOT USED, HAVING BEEN THE OBVIOUS CHOICE
%   -----------------------------------------------------
%   Levine 2008 is the famous number -- 53-57% myofiber CSA lost after 18-69
%   h of complete inactivity -- and calibrating k_deg from it was this
%   function's original design. It is wrong to do so. The 18-69 h window
%   spans a 3.8x range in implied k (0.28-1.07/day), which is 4-17x every
%   other study in the literature; it is a between-group biopsy comparison
%   at one time point, not a time course; the cases are brain-dead organ
%   donors carrying their own catabolic storm and the controls are thoracic
%   surgery patients ventilated 2-3 h. Levine establishes that rapid disuse
%   atrophy is real. It cannot fix a rate, and using it would have made
%   k_deg wrong by an order of magnitude.
%
%   WHAT REMAINS UNIDENTIFIED
%   -------------------------
%   g(1) -- the degradation at full spontaneous activity -- has no data
%   constraining it: Zambon's only measurement there is a net POSITIVE rate,
%   which the fit cannot use. The fit therefore drives g_min to its bound.
%   That is not a failure: g(1) = 0 says full activity leaves no net
%   degradation and capacity recovers to baseline, which is exactly what
%   Zambon's +2.3%/day at spontaneous breathing describes. It is reported
%   rather than hidden, and it is why k_syn must NOT be derived from g(1) --
%   an earlier version of this function did, and the circular dependency
%   collapsed k_syn to zero.
%
%   Name-value:
%       FitZambon     solve g_min/p_mono against Zambon (default true)
%       A_PS          activity assigned to 'pressure support' in Lecronier
%                     (default from config; OURS, not measured)
%
%   See also vidd.degradation, vidd.equilibriumCapacity, vidd.params

arguments
    cfg (1,1) struct
    opts.FitZambon (1,1) logical = true
    opts.A_PS      (1,1) double  = NaN
end

p = vidd.params(cfg);

% ---------------- g_min, p_mono: least squares on Zambon ----------------
z = cfg.validation.zambon2016;
A_z    = reshape(z.assigned_A, 1, []);
rate_z = reshape(z.dTdi_pct_per_day, 1, []) / 100;    % fractional /day

% Zambon's reported %/day is a fractional rate of change measured while
% patients are near their own baseline, so it reads as
%       (dC/dt)/C  ~  -(k_deg*g(A) + d) + k_syn*h(A)*(C_max0 - C)/C
% The three NEGATIVE points sit where synthesis is small (low A, h ~ h0), so
% they identify g up to a scale. The POSITIVE point at A=1 cannot: a positive
% net rate means synthesis dominates, and separating the two would need a
% measurement nobody has made. It is used for the k_syn solve instead.
%
% JABER SETS THE SCALE, ZAMBON SETS THE SHAPE -- and they must not be mixed.
% The two anchors disagree on the magnitude: Zambon's CMV rate is 0.078/day
% (thickness), Jaber's is 0.064/day (force), a 22% gap. Normalising Zambon's
% points by JABER's k_deg would force g(0) = 1.22 against a functional form
% that has g(0) = 1 by construction, and bury that disagreement in the fit
% residual -- which is what an earlier version of this function did, for an
% RMSE of 0.15 that was really a units conflict.
%
% Normalising Zambon's points by ZAMBON's OWN A=0 rate makes g(0) = 1 exactly
% and leaves the fit to do only what it can do: determine the SHAPE. The
% scale stays with Jaber, which is the better anchor (force, in the model's
% own units). The 22% disagreement is then visible in info.scale_conflict
% rather than hidden.
useForG  = rate_z < 0;
A_fit    = A_z(useForG);
rate_fit = rate_z(useForG);
g_target = rate_fit / rate_fit(A_fit == 0);   % normalised by Zambon's own CMV rate

info.scale_conflict_zambon_k = -rate_fit(A_fit == 0);
info.scale_conflict_jaber_k  = p.k_deg;
info.scale_conflict_ratio    = info.scale_conflict_zambon_k / p.k_deg;

if opts.FitZambon
    obj = @(x) sum((localGMono(A_fit, x(1), x(2)) - g_target).^2);
    x0  = [0.2, 1.5];
    lb  = [0, 0.2];
    ub  = [0.99, 8];
    optsFmin = optimoptions('fmincon', 'Display', 'none');
    xHat = fmincon(obj, x0, [], [], [], [], lb, ub, [], optsFmin);
    p.g_min  = xHat(1);
    p.p_mono = xHat(2);

    fitted = localGMono(A_fit, p.g_min, p.p_mono);
    info.zambon_A        = A_fit;
    info.zambon_g_target = g_target;
    info.zambon_g_fitted = fitted;
    info.zambon_rmse     = sqrt(mean((fitted - g_target).^2));
end

% ---------------- k_syn: Lecronier's two-group solve ----------------
L = cfg.validation.lecronier2022;
A_PS = opts.A_PS;
if isnan(A_PS), A_PS = L.assigned_A_at_second_measurement; end

C0s = L.septic_Ptr_stim_cmH2O;
C0n = L.nonseptic_Ptr_stim_cmH2O;
Rs  = 1 + L.septic_change_pct/100;
Rn  = 1 + L.nonseptic_change_pct/100;
tL  = L.median_days;

% Both groups share C* and tau at the same activity, so subtracting the two
% exponentials eliminates C*:
%   Rs*C0s = C* + (C0s - C*)*E ;  Rn*C0n = C* + (C0n - C*)*E
%   => (Rs*C0s - Rn*C0n) = E*(C0s - C0n)  ... after cancelling C*(1-E)
E = (Rs*C0s - Rn*C0n) / (C0s - C0n);
if ~(E > 0 && E < 1)
    error('vidd:calibrate:lecronierInfeasible', ...
        ['Lecronier two-point solve gives E = %.3f, outside (0,1): the two groups cannot be ' ...
         'described by one exponential relaxing to one C*.'], E);
end
tau_PS  = -tL / log(E);
Cstar_PS = (Rs*C0s - C0s*E) / (1 - E);

% tau = 1/(syn + loss) with loss = k_deg*g(A_PS) already known.
totalRate = 1 / tau_PS;
loss_PS   = p.k_deg * vidd.degradation(A_PS, p);
syn_PS    = totalRate - loss_PS;

if syn_PS <= 0
    error('vidd:calibrate:nonPositiveSynthesis', ...
        ['Lecronier''s relaxation rate (%.4f /day) is below the degradation rate at ' ...
         'A_PS = %.2f (%.4f /day), leaving no room for synthesis. Check the A assignment.'], ...
        totalRate, A_PS, loss_PS);
end

p.k_syn = syn_PS / vidd.synthesis(A_PS, p);

% The solve's own by-product: how far below baseline the septic patients
% started. C*/C_max0 = syn/(syn+loss), so C_max0 on Lecronier's scale is
% Cstar_PS/(syn/(syn+loss)).
C_max0_implied = Cstar_PS * totalRate / syn_PS;
info.lecronier_E              = E;
info.lecronier_tau_days       = tau_PS;
info.lecronier_Cstar          = Cstar_PS;
info.lecronier_C_max0_implied = C_max0_implied;
info.lecronier_septic_offset  = 1 - C0s / C_max0_implied;
info.demoule_septic_offset    = cfg.disease.sepsis_offset_fraction;
info.offset_agreement         = abs(info.lecronier_septic_offset - info.demoule_septic_offset);

% ---------------- report ----------------
info.k_deg  = p.k_deg;
info.k_syn  = p.k_syn;
info.g_min  = p.g_min;
info.p_mono = p.p_mono;
info.g_mode = p.g_mode;
info.d_mode = p.d_mode;
info.A_PS   = A_PS;

info.C_star_at_1_nodisease = vidd.equilibriumCapacity(1, setD(p, 0));
info.C_star_at_0_nodisease = vidd.equilibriumCapacity(0, setD(p, 0));
info.k_syn_identified      = true;    % from Lecronier's two groups

% What Jaber should reproduce: -32% force at 6 days under ventilation.
info.jaber = cfg.validation.jaber2010;
info.jaber_predicted_loss_at_A0 = 1 - localFracAt(p, 0, cfg.validation.jaber2010.duration_days);

info.constraints = [
    "k_deg  <- ANCHORED: Jaber 2010 PMID 20813887, TwPtr -32% at 6 d -> 0.064/day"
    "g_min, p_mono <- SOLVED: least squares on Zambon 2016 PMID 26992064 (negative-rate points only)"
    "k_syn  <- SOLVED: Lecronier 2022 PMID 35403916, two groups at one intervention"
];

info.notes = [
    "Levine 2008 is deliberately EXCLUDED from the rate calibration: its 18-69 h window implies k = 0.28-1.07/day, 4-17x every other study."
    "g(1) is unconstrained by data -- Zambon's only measurement at full activity is a net POSITIVE rate the fit cannot use. g_min lands on its bound, meaning full activity leaves no net degradation, which is what Zambon describes."
    "Zambon's mode -> A assignment and Lecronier's A_PS are OURS, not measured. Where the fitted g sits on the A axis inherits that."
    "Capacity is read as proportional to cross-sectional area, and rates measured on twitch pressure are transferred to MIP as fractional changes. Both are declared assumptions."
];
end

% ---------------- helpers ----------------

function g = localGMono(A, g_min, p_mono)
g = g_min + (1 - g_min) .* (1 - A) .^ p_mono;
end

function frac = localFracAt(p, A, tDays)
% C(t)/C_max0 starting from the ceiling, for the linear system.
syn  = p.k_syn * vidd.synthesis(A, p);
loss = p.k_deg * vidd.degradation(A, p) + p.d_disease;
r = syn / (syn + loss);
S = syn + loss;
frac = r + (1 - r) * exp(-S * tDays);
end

function q = setD(p, d)
q = p; q.d_disease = d;
end
