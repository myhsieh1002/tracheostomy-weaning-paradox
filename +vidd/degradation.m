function g = degradation(A, p)
%DEGRADATION Dimensionless degradation multiplier vs diaphragm activity.
%
%   g = vidd.degradation(A, p) returns the multiplier on the baseline
%   degradation rate k_deg. By convention g(0) = 1, so k_deg IS the
%   complete-disuse rate and carries the units.
%
%   TWO MODES, BECAUSE THE LITERATURE SUPPORTS TWO DIFFERENT SHAPES
%   --------------------------------------------------------------
%   p.g_mode = 'monotonic'  (default -- what the rate data show)
%
%       g(A) = g_min + (1 - g_min)*(1 - A)^p_mono
%
%       Calibrated against Zambon 2016 (PMID 26992064), the ONLY study
%       stratifying diaphragm atrophy RATE by ventilatory support:
%       -7.5%/d (CMV), -5.3% (high PSV), -1.5% (low PSV), +2.3%
%       (spontaneous/CPAP). Monotonic, best at the fully spontaneous end,
%       no upturn. Zambon's own conclusion: "a linear relationship between
%       ventilator support and diaphragmatic atrophy rate."
%
%   p.g_mode = 'ushape'  (the spec's version)
%
%       g(A) = g0 + a_disuse*(A* - A)_+^p + a_injury*(A - A*)_+^q
%
%   WHY MONOTONIC IS THE DEFAULT, WHEN THE U IS THE SPEC'S HYPOTHESIS
%   ----------------------------------------------------------------
%   Because the U is measured at a different level than this function
%   operates at. Goligher 2018 (PMID 28930478, n=191) is real and
%   well-adjusted: a thickening fraction of 15-30% gives the shortest time
%   to liberation, both tails worse. But that is a U in a CLINICAL OUTCOME.
%   g is a rate of capacity loss, and nobody has measured a U in that.
%   The only study that looked (Zambon) found the opposite.
%
%   The two are reconcilable -- if high-effort injury manifests as
%   dysfunctional THICKENING rather than thinning, atrophy rate could be
%   monotonic while function is U-shaped. That is Goligher's own reading and
%   is consistent with maximal thickening fraction being depressed in BOTH
%   tails (Goligher 2015). It has never been quantified.
%
%   So: 'monotonic' is what the rate data support; 'ushape' is the
%   hypothesis. Running both, and reporting which conclusions survive, is
%   the honest treatment -- see experiments/vidd_exp1_shape.m.
%
%   THE SHAPE OF THE U IS NOT DATA
%   ------------------------------
%   In 'ushape' mode, a_injury has NO published constraint at all:
%   Orozco-Levi 2001 proves load-induced sarcomere disruption exists in
%   humans but gives no dose-response, no time course and no force
%   measurement. The exponents p and q are pure invention -- no study has
%   measured the shape of degradation against effort. They are sensitivity
%   parameters and are reported as such.
%
%   See also vidd.synthesis, vidd.calibrate, vidd.dCdt

arguments
    A double
    p (1,1) struct
end

if any(A < 0 | A > 1, 'all')
    error('vidd:degradation:activityOutOfRange', ...
        'A must lie in [0, 1]; got values in [%g, %g].', min(A(:)), max(A(:)));
end

switch lower(string(p.g_mode))
    case "monotonic"
        g = p.g_min + (1 - p.g_min) .* (1 - A) .^ p.p_mono;

    case "ushape"
        disuse = max(p.A_star - A, 0) .^ p.p_exp;
        injury = max(A - p.A_star, 0) .^ p.q_exp;
        g = p.g0 + p.a_disuse .* disuse + p.a_injury .* injury;

    otherwise
        error('vidd:degradation:unknownMode', ...
            'g_mode must be ''monotonic'' or ''ushape'' (got "%s").', p.g_mode);
end
end
