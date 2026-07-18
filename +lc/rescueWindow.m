function rw = rescueWindow(C_max, L_device_range, p)
%RESCUEWINDOW Named integration API -- fold structure and rescue window at one C_max.
%
%   rw = lc.rescueWindow(C_max, L_device_range, p) is the stable, named
%   interface specified in Model 2's build plan section 2.5b. It is the
%   single seam through which Model 2b (VIDD capacity dynamics) drives this
%   model: 2b supplies a capacity trajectory C_max(t) and calls this
%   function once per time point to recompute the folds. Nothing else may
%   re-derive the fold; there must be exactly one implementation of it.
%
%   INPUTS
%       C_max           scalar capacity ceiling (cmH2O). For the dynamic
%                       mode, pass C_max(t) one point at a time -- see
%                       lc.rescueWindowTrajectory.
%       L_device_range  [L_low, L_high] the achievable load span, i.e.
%                       [L_TRACH, L_ETT] from Model 1 via lc.coupling.
%                       NEVER hard-code these.
%       p               dynamical parameters (lc.params / lc.calibrate)
%
%   OUTPUT struct `rw`
%       fold_left    saddle-node at lower load (cmH2O); recovery boundary
%       fold_right   saddle-node at higher load (cmH2O); collapse boundary
%       bistable     true when both folds exist at physical load
%       window       [L_lo, L_hi]: the sub-interval of L_device_range that
%                    lands on the sustainable branch. Empty when no
%                    achievable load sustains the patient.
%       spans_fold   true when L_device_range straddles fold_right -- i.e.
%                    the device choice DECIDES the outcome. This is the
%                    rescue condition.
%       basin_at     function handle (C, L) -> "high" | "low", giving the
%                    attractor a state converges to
%       C_max_window the analytic band of C_max over which this
%                    L_device_range yields a rescue: (L_low/l_high,
%                    L_high/l_high)
%
%   THE UNDERLYING CLOSED FORM
%   --------------------------
%   By the scale invariance (lc.normalizedFolds), folds sit at fixed
%   normalised load and scale linearly:  L_fold = l_fold * C_max. So this
%   function does no root-finding per call beyond the one-off normalised
%   fold computation, which makes the per-time-point dynamic mode cheap.
%
%   See also lc.normalizedFolds, lc.rescueOutcome, lc.rescueWindowTrajectory

arguments
    C_max          (1,1) double {mustBePositive}
    L_device_range (1,2) double {mustBeNonnegative}
    p              (1,1) struct
end

if L_device_range(1) > L_device_range(2)
    error('lc:rescueWindow:badRange', ...
        'L_device_range must be [low, high]; got [%g, %g].', L_device_range(1), L_device_range(2));
end

nf = lc.normalizedFolds(p);

rw.bistable = nf.isBistable;

if ~nf.isBistable
    % No fold: the system is monostable and the notion of a rescue does not
    % apply. Return a well-formed struct rather than erroring, so a caller
    % sweeping C_max(t) can record "no window" without special-casing.
    rw.fold_left    = NaN;
    rw.fold_right   = NaN;
    rw.window       = [];
    rw.spans_fold   = false;
    rw.C_max_window = [NaN NaN];
    rw.basin_at     = @(C, L) basinAt(C, L, C_max, p);
    rw.C_max        = C_max;
    rw.l_high       = NaN;
    rw.l_low        = NaN;
    return;
end

rw.fold_left  = nf.l_low  * C_max;
rw.fold_right = nf.l_high * C_max;
rw.l_low      = nf.l_low;
rw.l_high     = nf.l_high;
rw.C_max      = C_max;

% Loads below fold_right leave a sustainable branch in existence.
lo = L_device_range(1);
hi = min(L_device_range(2), rw.fold_right);
if lo < rw.fold_right
    rw.window = [lo, hi];
else
    rw.window = [];   % even the lightest achievable load is above the fold
end

% The device decides the outcome only when the two ends straddle the fold.
rw.spans_fold = L_device_range(1) < rw.fold_right && L_device_range(2) > rw.fold_right;

% The C_max band over which this load pair produces a rescue.
rw.C_max_window = [L_device_range(1) / nf.l_high, L_device_range(2) / nf.l_high];

rw.basin_at = @(C, L) basinAt(C, L, C_max, p);
end

function which = basinAt(C, L, C_max, p)
%BASINAT Which attractor does state C converge to under load L?
%
%   The separatrix is the unstable middle equilibrium: above it the
%   trajectory rises to the high branch, below it falls to the low one.
%   Resolving it by the fixed-point structure rather than by integrating is
%   exact and has no horizon to choose.

fp = lc.fixedPoints(L, C_max, p);

if numel(fp.C) == 1
    % Monostable: everything goes to the sole equilibrium.
    which = ternary(fp.C(1) > 0.5 * C_max, "high", "low");
    return;
end

separatrix = fp.C(2);   % the middle root is the unstable one
which = ternary(C > separatrix, "high", "low");
end

function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end
