function [p, info] = calibrate(cfg, opts)
%CALIBRATE Pin down beta and s against published physiology.
%
%   [p, info] = lc.calibrate(cfg) returns a parameter struct in which the
%   two parameters that have no directly published value -- beta and s --
%   are determined by two independent literature constraints, rather than
%   left as placeholders.
%
%   THE TWO CONSTRAINTS
%   -------------------
%   The model has four parameters. Two are anchored directly:
%
%       alpha  = 0.15 /hr   Laghi 1995 (PMID 7592215), 24 h recovery curve
%       u_crit = 0.4        Roussos & Macklem 1977 (PMID 893274), Pdicrit
%
%   The remaining two are fixed by two measured quantities that the model
%   reproduces analytically:
%
%   (1) DEPTH OF THE FAILURE ATTRACTOR  ->  beta
%       The low branch satisfies, exactly,
%
%           C_low / C_max = alpha / (alpha + beta)
%
%       Laghi 1995 fatigued healthy subjects to task failure and measured
%       twitch Pdi falling 38.9 -> 25.1 cmH2O, i.e. to 64.5% of baseline.
%       Setting C_low/C_max = 0.645 gives
%
%           beta = alpha * (1 - 0.645) / 0.645
%
%   (2) LOCATION OF THE RIGHT FOLD  ->  s
%       By the scale invariance (lc.normalizedFolds), the fold sits at a
%       fixed normalised load l_high = L/C_max, which IS Pi/Pimax -- given
%       that Pimax is the rested value measured before the trial, so that it
%       reads as the capacity ceiling rather than an already-fatigued state.
%       Vassilakopoulos 1998 (PMID 9700110) measured, in the same patients,
%       Pi/Pimax = 0.31 when they weaned and 0.46 when they failed. For the
%       model to reproduce that discrimination, the fold must separate them:
%
%           0.310  <  l_high  <  0.461
%
%       s is solved so that l_high hits the target.
%
%   HOW WELL IS s ACTUALLY IDENTIFIED? (weakly -- state this in the paper)
%   ---------------------------------------------------------------------
%   l_high is bounded above by u_crit and approaches it only as s -> inf, so
%   the reachable part of the discrimination band is (0.310, 0.400), not the
%   full (0.310, 0.461). Every s above ~23 lands inside it:
%
%       l_high = 0.310  ->  s ~ 23
%       l_high = 0.350  ->  s ~ 71
%       l_high = 0.385  ->  s ~ 354
%       l_high -> 0.400 ->  s -> inf
%
%   So the discrimination constrains s only from BELOW. Near the ceiling a
%   large change in s barely moves the fold, which is why targeting the
%   midpoint of the Vassilakopoulos band (0.385) forces an implausibly hard
%   threshold (s ~ 354, transition width du ~ 0.011) -- that is
%   over-reading the data. The default target is instead the midpoint of the
%   REACHABLE band, 0.355, giving s ~ 76. s remains a sensitivity parameter.
%
%   WHY THIS MATTERS
%   ----------------
%   The build spec carried beta = 0.8 and s = 25 as placeholders. With
%   those values l_high = 0.218, and BOTH Vassilakopoulos groups fall above
%   the fold -- the model predicts failure for the patients who actually
%   weaned. The calibrated values reproduce the observed discrimination.
%
%   A caveat that belongs in the manuscript: constraint (1) comes from
%   healthy subjects in an experimental fatigue protocol, and equating that
%   capacity drop with the model's low attractor is an interpretation, not a
%   measurement. It is nonetheless the only quantitative anchor available
%   for beta, and it is far better than an arbitrary value.
%
%   A SECOND, LOAD-BEARING INVARIANCE
%   ---------------------------------
%   The equilibrium condition reduces to g(u) = (alpha/beta)*(1-x)/x, which
%   contains alpha and beta ONLY as their ratio. The entire bifurcation
%   structure therefore depends on (beta/alpha, s, u_crit); the absolute
%   rates set the timescale of approach and nothing else. Verified over a
%   100-fold range of alpha in tests/tLoadCapacity.m.
%
%   This retires the calibration's biggest apparent weakness. Laghi's
%   recovery curve is biphasic and no single alpha fits it -- but alpha's
%   absolute value does not enter H1, H2 or H3 at all. It matters only for
%   M2b (time-to-failure and the early-warning signals), where the timescale
%   is the quantity of interest.
%
%   Name-value options:
%       TargetLHigh   normalised fold to solve s for
%                     (default: midpoint of the REACHABLE band, 0.355)
%       CapacityDrop  C_low/C_max from Laghi (default 0.645)
%       SBracket      search bracket for s (default [1, 5000])
%
%   `info` records the constraints, the solved values, the achieved
%   l_high/l_low, and whether the discrimination check passes.
%
%   See also lc.normalizedFolds, lc.params

arguments
    cfg (1,1) struct
    opts.TargetLHigh  (1,1) double = NaN
    opts.CapacityDrop (1,1) double {mustBeInRange(opts.CapacityDrop, 0, 1)} = 0.645
    opts.SBracket     (1,2) double = [1, 5000]
end

v = cfg.validation.vassilakopoulos1998;
l_success = v.success.u;
l_failure = v.failure.u;

p = lc.params(cfg);

if isnan(opts.TargetLHigh)
    % Midpoint of the REACHABLE band: the discrimination requires
    % l_high > l_success, and the structure caps it below u_crit.
    targetLHigh = mean([l_success, p.u_crit]);
else
    targetLHigh = opts.TargetLHigh;
end

% --- Constraint (1): beta from the failure-attractor depth ---
frac = opts.CapacityDrop;
p.beta = p.alpha * (1 - frac) / frac;

% --- Constraint (2): s from the fold location ---
objective = @(s) foldResidual(s, p, targetLHigh);

fa = objective(opts.SBracket(1));
fb = objective(opts.SBracket(2));
if ~isfinite(fa) || ~isfinite(fb) || fa * fb > 0
    error('lc:calibrate:noBracket', ...
        ['Cannot bracket s in [%g, %g] for target l_high = %.4f. ' ...
         'l_high is bounded above by u_crit = %.2f, so no s can reach a target at or above it.'], ...
        opts.SBracket(1), opts.SBracket(2), targetLHigh, p.u_crit);
end

p.s = fzero(objective, opts.SBracket);

nf = lc.normalizedFolds(p);

info.alpha        = p.alpha;
info.beta         = p.beta;
info.s            = p.s;
info.u_crit       = p.u_crit;
info.targetLHigh  = targetLHigh;
info.achievedLHigh = nf.l_high;
info.achievedLLow  = nf.l_low;
info.C_low_frac    = nf.C_low_frac;
info.capacityDrop  = frac;
info.l_success     = l_success;
info.l_failure     = l_failure;

% The discrimination the calibration exists to reproduce.
info.discriminates = nf.l_high > l_success && nf.l_high < l_failure;
info.predictSuccessWeans = l_success < nf.l_high;
info.predictFailureFails = l_failure > nf.l_high;

info.constraints = [
    "beta  <- Laghi 1995 PMID 7592215: C_low/C_max = " + string(frac)
    "s     <- Vassilakopoulos 1998 PMID 9700110: l_high = " + string(targetLHigh)
    "alpha <- Laghi 1995 PMID 7592215 (24h recovery)"
    "u_crit<- Roussos & Macklem 1977 PMID 893274 (Pdicrit)"
];
end

function r = foldResidual(s, p, target)
p.s = s;
nf = lc.normalizedFolds(p);
if ~nf.isBistable
    r = -target;   % no fold: drive the search towards larger s
else
    r = nf.l_high - target;
end
end
