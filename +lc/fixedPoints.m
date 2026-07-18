function fp = fixedPoints(L, C_max, p, opts)
%FIXEDPOINTS Locate every equilibrium of the load-capacity system.
%
%   fp = lc.fixedPoints(L, C_max, p) returns a struct with the equilibria
%   and their stability for a given load L and capacity ceiling C_max.
%
%   fp = lc.fixedPoints(..., ScanN=n, CMinFrac=r) controls the bracketing
%   scan resolution and its lower bound (as a fraction of C_max).
%
%   METHOD
%   ------
%   All equilibria lie in the open interval (0, C_max):
%
%     * as C -> 0+,  f -> alpha*C_max > 0
%     * at  C = C_max, f = -beta*C_max*g(L/C_max) < 0
%     * for C > C_max, both terms are negative, so f < 0
%
%   f is continuous there, so sign changes bracket every root and their
%   number is necessarily ODD -- one or three. A dense scan brackets them
%   and fzero refines each to machine precision.
%
%   A bracketing scan is used in preference to multi-start fsolve because
%   it cannot miss a root that lies between two scan points *and* report
%   success: the parity check below turns any missed pair into a hard error
%   rather than a silently monostable answer. Near a fold the two outer
%   roots approach each other, which is exactly where a Newton method
%   started from a fixed guess set tends to converge onto the same root
%   twice and quietly report bistability as monostability.
%
%   Fields of `fp`:
%       C          equilibria, ascending (1x1 or 1x3)
%       stable     logical, true where df/dC < 0
%       lambda     df/dC at each equilibrium
%       nStable    number of stable equilibria
%       isBistable true when there are two stable branches
%       branch     "low" | "middle" | "high" label per equilibrium
%       L, C_max   echoed inputs
%
%   See also lc.dCdt, lc.continuation, lc.foldPoints

arguments
    L     (1,1) double {mustBeNonnegative}
    C_max (1,1) double {mustBePositive}
    p     (1,1) struct
    opts.ScanN    (1,1) double {mustBePositive} = 2000
    opts.CMinFrac (1,1) double {mustBePositive} = 1e-4
    opts.RootTol  (1,1) double {mustBePositive} = 1e-12
end

Cgrid = linspace(C_max * opts.CMinFrac, C_max, opts.ScanN);
fvals = lc.dCdt(Cgrid, L, C_max, p);

% Bracket every sign change.
crossings = find(fvals(1:end-1) .* fvals(2:end) < 0);

roots = zeros(1, numel(crossings));
fzOpts = optimset('TolX', opts.RootTol);
for k = 1:numel(crossings)
    i = crossings(k);
    roots(k) = fzero(@(C) lc.dCdt(C, L, C_max, p), [Cgrid(i), Cgrid(i+1)], fzOpts);
end

% Catch an exact zero landing on a scan point, which produces no sign change.
exactZeros = Cgrid(fvals == 0);
roots = unique([roots, exactZeros]);

if isempty(roots)
    error('lc:fixedPoints:noRootFound', ...
        ['No equilibrium found for L=%g, C_max=%g. The vector field must change ' ...
         'sign on (0, C_max]; this indicates a parameter or scan-bound problem.'], ...
        L, C_max);
end

% The root count must be odd. An even count means the scan stepped over a
% pair of roots -- fail loudly rather than report a wrong bifurcation
% structure.
if mod(numel(roots), 2) == 0
    error('lc:fixedPoints:evenRootCount', ...
        ['Found %d equilibria for L=%g, C_max=%g; the count must be odd. ' ...
         'The scan (ScanN=%d) has stepped over a closely spaced pair -- increase ScanN.'], ...
        numel(roots), L, C_max, opts.ScanN);
end

[~, dfdC] = lc.dCdt(roots, L, C_max, p);

fp.C          = roots;
fp.lambda     = dfdC;
fp.stable     = dfdC < 0;
fp.nStable    = sum(fp.stable);
fp.isBistable = numel(roots) == 3 && fp.nStable == 2;
fp.L          = L;
fp.C_max      = C_max;

fp.branch = strings(1, numel(roots));
if numel(roots) == 1
    fp.branch(1) = "single";
elseif numel(roots) == 3
    fp.branch = ["low", "middle", "high"];
else
    fp.branch(:) = "unlabelled";
end
end
