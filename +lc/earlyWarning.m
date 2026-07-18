function ew = earlyWarning(t, C, opts)
%EARLYWARNING Rolling variance and lag-1 autocorrelation as fold indicators.
%
%   ew = lc.earlyWarning(t, C) computes the two classic critical-slowing-down
%   statistics over a sliding window.
%
%   Name-value options:
%       Window     window length in samples of the ANALYSED (subsampled)
%                  series (default 10% of it)
%       Detrend    remove a rolling mean before computing statistics
%                  (default true)
%       Subsample  keep every n-th sample before analysing (default 1)
%
%   WHY SUBSAMPLE MATTERS FOR AR(1)
%   -------------------------------
%   AR(1) measures how much of a perturbation survives ONE sampling
%   interval. Integrating with a step far shorter than the system's
%   relaxation time makes that answer trivially "almost all of it": with
%   dt = 0.01 h against a relaxation time of ~1/alpha ~ 7 h, AR(1) pins at
%   ~0.998 everywhere and has no dynamic range left to report the slowing
%   down. The integrator needs the small step for accuracy; the STATISTIC
%   needs a sampling interval comparable to the relaxation time. Subsample
%   decouples the two -- pick it so that Subsample*dt is a useful fraction
%   of 1/alpha.
%
%   WHY DETRENDING IS NOT OPTIONAL IN PRACTICE
%   ------------------------------------------
%   As the load ramps, C drifts downward along the branch. That drift is a
%   deterministic trend, not a fluctuation, and it inflates the rolling
%   variance on its own -- so an undetrended "variance rises near the fold"
%   result would be partly an artefact of the ramp rather than evidence of
%   critical slowing down. Subtracting a rolling mean first leaves the
%   residual fluctuations, which is what the theory is about.
%
%   AR(1) is estimated on the same residuals by ordinary least squares
%   (regressing x(k+1) on x(k)); as a fold is approached the restoring
%   eigenvalue goes to zero, perturbations decay more slowly, and AR(1)
%   rises towards 1.
%
%   Fields of `ew`: tMid, variance, ar1, window, residual.
%
%   See also lc.simulateStochastic

arguments
    t (1,:) double
    C (1,:) double
    opts.Window    (1,1) double = NaN
    opts.Detrend   (1,1) logical = true
    opts.Subsample (1,1) double {mustBePositive} = 1
end

if opts.Subsample > 1
    keep = 1:round(opts.Subsample):numel(C);
    t = t(keep);
    C = C(keep);
end

n = numel(C);
w = opts.Window;
if isnan(w), w = max(20, round(0.10*n)); end
w = min(w, n);

if opts.Detrend
    % Rolling mean over the same window; 'shrink' keeps the endpoints
    % defined rather than returning NaN there.
    trend = movmean(C, w, 'Endpoints','shrink');
    resid = C - trend;
else
    resid = C - mean(C);
end

nWin = n - w + 1;
tMid = zeros(1, nWin);
variance = zeros(1, nWin);
ar1 = zeros(1, nWin);

for k = 1:nWin
    idx = k:(k+w-1);
    seg = resid(idx);
    tMid(k) = t(idx(round(w/2)));
    variance(k) = var(seg);

    x = seg(1:end-1)';
    y = seg(2:end)';
    x = x - mean(x); y = y - mean(y);
    denom = sum(x.^2);
    if denom > 0
        ar1(k) = sum(x .* y) / denom;
    else
        ar1(k) = NaN;
    end
end

ew.tMid     = tMid;
ew.variance = variance;
ew.ar1      = ar1;
ew.window   = w;
ew.residual = resid;
end
