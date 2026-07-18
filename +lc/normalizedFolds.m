function nf = normalizedFolds(p, opts)
%NORMALIZEDFOLDS Scale-invariant fold locations of the load-capacity system.
%
%   nf = lc.normalizedFolds(p) returns the saddle-node folds in NORMALISED
%   load, l = L/C_max, which is the only place they live.
%
%   THE SCALE INVARIANCE
%   --------------------
%   Substituting C = C_max*x and L = C_max*l into the vector field,
%
%       dC/dt = alpha*(C_max - C) - beta*C*g(L/C)
%     =>  dx/dt = alpha*(1 - x) - beta*x*g(l/x)
%
%   C_max cancels EXACTLY. The dynamics depend only on the normalised
%   capacity x and the normalised load l -- never on C_max itself. So:
%
%     * the folds sit at fixed l, and move linearly in absolute load:
%       L_fold = l_fold * C_max;
%     * l = L/C_max is precisely Pi/Pimax, the ratio the weaning
%       literature reports, so the model's fold is directly comparable to
%       published thresholds;
%     * the bistable window in absolute load NARROWS in proportion to
%       C_max, which is the analytic core of the rescue-window result.
%
%   This is verified numerically over a 200-fold range of C_max in
%   tests/tLoadCapacity.m.
%
%   RELATION TO u_crit
%   ------------------
%   l_high < u_crit ALWAYS, approaching it from below as s -> inf. The fold
%   is not at u_crit: it is where beta*C*g(u) first overtakes
%   alpha*(C_max - C), which happens while g is still well below 0.5. The
%   spec's sanity check -- that a hard threshold drives the fold to u_crit --
%   holds, and is asserted in the tests.
%
%   Fields of `nf`: l_low, l_high, x_low, x_high, isBistable, width,
%   C_low_frac (= alpha/(alpha+beta), the depth of the failure attractor).
%
%   See also lc.equilibriumBranch, lc.calibrate

%   PERFORMANCE
%   -----------
%   This is the innermost function of every Model 2 analysis, and the Sobol
%   designs call it once per sample with a DIFFERENT parameter vector, so it
%   cannot be hoisted out of those loops the way m2_exp3 hoists it. It
%   therefore solves the folds directly -- a coarse bracketing scan of
%   dl/dx followed by fzero -- instead of building the full 20,000-point
%   branch via lc.equilibriumBranch and reading the folds off it. Same
%   answer to fzero's tolerance, ~100x faster.
%
%   tests/tLoadCapacity.m asserts this agrees with the branch-derived folds.

arguments
    p (1,1) struct
    opts.ScanN (1,1) double {mustBePositive} = 400
end

% C_max is arbitrary by the invariance above; 1 makes L numerically equal
% to l, and x = C.
xlo = p.alpha / (p.alpha + p.beta);
xhi = 1;

% As beta -> 0 the domain (alpha/(alpha+beta), 1) collapses towards a point.
% Once its width approaches double-precision resolution near 1, the inset
% used to step inside it is smaller than eps(1) and G = alpha*(1-x)/(beta*x)
% rounds to values outside (0,1) -- making log(G/(1-G)) complex. That is a
% floating-point failure of the DOMAIN, not a fold. Bail out early: no
% fatigue means no fold, which is exactly the beta -> 0 limit the spec asks
% for.
if (xhi - xlo) < 1e-6
    nf = noFold(p);
    return;
end

eps_ = max((xhi - xlo) * 1e-9, eps(1));
x = linspace(xlo + eps_, xhi - eps_, opts.ScanN);

% Keep only samples on which the branch is genuinely defined.
G = p.alpha * (1 - x) ./ (p.beta * x);
valid = isfinite(G) & G > 0 & G < 1;
x = x(valid);
if numel(x) < 3
    nf = noFold(p);
    return;
end

d = dldx(x, p);

% Folds are the stationary points of l(x): sign changes of dl/dx. Bracket
% only across pairs where both endpoints are finite.
finiteOK = isfinite(d);
idx = find(d(1:end-1) .* d(2:end) < 0 & finiteOK(1:end-1) & finiteOK(2:end));

folds = zeros(1, numel(idx));
for k = 1:numel(idx)
    folds(k) = fzero(@(xx) dldx(xx, p), [x(idx(k)), x(idx(k)+1)]);
end

if numel(folds) == 2
    lvals = lOf(folds, p);
    [nf.l_low, iLo] = min(lvals);
    [nf.l_high, iHi] = max(lvals);
    nf.isBistable = nf.l_high > 0;
    nf.width      = nf.l_high - nf.l_low;
    nf.x_low      = folds(iLo);
    nf.x_high     = folds(iHi);
else
    % Below the cusp the two folds have merged: no bistability.
    nf.isBistable = false;
    nf.l_low  = NaN;
    nf.l_high = NaN;
    nf.width  = 0;
    nf.x_low  = NaN;
    nf.x_high = NaN;
end

% Depth of the low (failure) attractor. As C falls, u = L/C rises and
% g -> 1, so the balance alpha*(1-x) = beta*x gives x -> alpha/(alpha+beta)
% independently of load. This is an exact limit, and it is what lets Laghi's
% measured capacity drop pin down beta.
nf.C_low_frac = p.alpha / (p.alpha + p.beta);
end

function nf = noFold(p)
%NOFOLD Well-formed "monostable" result, for callers sweeping parameters.
nf.isBistable = false;
nf.l_low      = NaN;
nf.l_high     = NaN;
nf.width      = 0;
nf.x_low      = NaN;
nf.x_high     = NaN;
nf.C_low_frac = p.alpha / (p.alpha + p.beta);
end

% --- normalised branch: l(x) and dl/dx, with C_max = 1 ---

function l = lOf(x, p)
G = p.alpha * (1 - x) ./ (p.beta * x);
l = x .* (p.u_crit + log(G ./ (1 - G)) / p.s);
end

function d = dldx(x, p)
G = p.alpha * (1 - x) ./ (p.beta * x);
phi = log(G ./ (1 - G));
d = p.u_crit + phi / p.s - p.alpha ./ (p.s * p.beta * x .* G .* (1 - G));
end
