function lev = strategyLeverage(p, opts)
%STRATEGYLEVERAGE How much can ventilation strategy still buy? (H3)
%
%   lev = vidd.strategyLeverage(p) sweeps activity over [0,1] and returns
%   the best and worst equilibrium capacity a strategy can reach, and the
%   gap between them.
%
%   THIS IS H3, MADE QUANTITATIVE
%   ----------------------------
%   d_disease is independent of A, so it sets a floor on the loss rate that
%   no strategy can lift:
%
%       C*(A) = k_syn*h(A) / (k_deg*g(A) + d_disease)
%
%   The best a strategy can do is drive g(A) to its floor. As d_disease
%   grows it dominates the denominator, every A gives nearly the same
%   answer, and the ABSOLUTE leverage -- C*_best - C*_worst, in cmH2O --
%   collapses towards zero. Optimising activity stops buying anything.
%
%   ABSOLUTE AND RELATIVE LEVERAGE DISAGREE, AND THE DIFFERENCE IS REAL
%   ------------------------------------------------------------------
%   The RATIO C*_best/C*_worst does NOT go to 1 as d_disease grows -- it
%   goes to 1/h0, because the training term h(A) is a multiplicative factor
%   that survives in the numerator however large the denominator gets. Only
%   with h0 = 1 (h constant) does the ratio also collapse.
%
%   So the honest statement of H3 is about the ABSOLUTE gap: a catabolic
%   patient's capacity becomes insensitive to strategy in cmH2O, which is
%   the unit that matters, because it is the unit Model 2's fold lives in. A
%   30% relative advantage over a capacity of 8 cmH2O is not a rescue.
%   `lev.absolute` is the H3 metric; `lev.relative` is reported alongside so
%   the distinction is visible rather than buried in a choice of axis.
%
%   Name-value:
%       N          activity grid resolution (default 2001)
%       Threshold  a capacity of interest (e.g. where Model 2's rescue
%                  window sits); lev.bestAboveThreshold reports whether the
%                  best achievable C* clears it.
%
%   Fields of `lev`: A, Cstar, A_best, C_best, A_worst, C_worst, absolute,
%   relative, bestAboveThreshold, threshold.
%
%   See also vidd.equilibriumCapacity, vidd.dCdt

arguments
    p (1,1) struct
    opts.N         (1,1) double {mustBePositive} = 2001
    opts.Threshold (1,1) double = NaN
end

A = linspace(0, 1, opts.N);
Cstar = vidd.equilibriumCapacity(A, p);

[C_best,  iB] = max(Cstar);
[C_worst, iW] = min(Cstar);

lev.A       = A;
lev.Cstar   = Cstar;
lev.A_best  = A(iB);
lev.C_best  = C_best;
lev.A_worst = A(iW);
lev.C_worst = C_worst;

lev.absolute = C_best - C_worst;                       % cmH2O -- the H3 metric
lev.relative = C_best / max(C_worst, eps);             % dimensionless

lev.threshold = opts.Threshold;
if isnan(opts.Threshold)
    lev.bestAboveThreshold = NaN;
else
    lev.bestAboveThreshold = C_best > opts.Threshold;
end
end
