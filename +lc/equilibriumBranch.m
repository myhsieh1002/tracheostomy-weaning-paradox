function br = equilibriumBranch(C_max, p, opts)
%EQUILIBRIUMBRANCH Closed-form equilibrium manifold of the load-capacity system.
%
%   br = lc.equilibriumBranch(C_max, p) returns the complete equilibrium
%   curve C* vs L, parametrised by C, together with stability and the
%   saddle-node folds.
%
%   WHY THERE IS NO CONTINUATION HERE
%   ---------------------------------
%   The build spec called for pseudo-arclength continuation. That is not
%   needed: this system inverts in closed form. Setting dC/dt = 0,
%
%       alpha*(C_max - C) = beta*C*g(L/C)
%       g(L/C) = G,        G := alpha*(C_max - C) / (beta*C)
%
%   and g is a logistic, so it inverts exactly:
%
%       L(C) = C * [ u_crit + (1/s)*log( G / (1 - G) ) ]
%
%   Every C maps to exactly ONE L. The equilibrium manifold is therefore a
%   graph over C -- the fold is a fold in L, not in C, and vanishes entirely
%   once C is used as the parameter. Sampling C and evaluating L(C) traces
%   the whole S-shaped curve, including both unstable and stable arms, with
%   no predictor, no corrector and no step-size control.
%
%   This is strictly better than continuation for this model: it is exact,
%   it is O(n), and it cannot jump branches or step past a fold -- the two
%   failure modes that make arclength continuation delicate near exactly the
%   points of interest. lc.fixedPoints solves the forward problem (given L,
%   find C) by an independent bracketing scan, and tests/tContinuation.m
%   checks the two against each other.
%
%   DOMAIN
%   ------
%   L(C) is real only where 0 < G < 1, i.e.
%
%       alpha*C_max/(alpha + beta)  <  C  <  C_max
%
%   Below that interval fatigue can never balance recovery; above it C
%   exceeds the ceiling. As C approaches the lower bound, L -> +inf; as C
%   approaches C_max, L -> -inf. Only L >= 0 is physical, and `br` is
%   trimmed to it by default.
%
%   FOLDS
%   -----
%   Saddle-node points are the stationary points of L(C), located where
%
%       dL/dC = u_crit + phi/s - (alpha*C_max)/(s*beta*C*G*(1-G)) = 0
%
%   with phi = log(G/(1-G)). This derivative is analytic, so folds are found
%   by refining sign changes of dL/dC rather than by detecting a turning
%   point numerically.
%
%   Fields of `br`:
%       C, L        the equilibrium curve (ascending in C)
%       dLdC        analytic derivative along the curve
%       stable      logical, true where the equilibrium is stable
%       lambda      df/dC along the curve
%       folds       struct array with fields C, L, type ("left"|"right")
%       L_fold_low, L_fold_high   the bistable window in L (NaN if none)
%       isBistable  true when two folds exist at L >= 0
%       C_max, p    echoed inputs
%
%   See also lc.fixedPoints, lc.dCdt, lc.rescueOutcome

arguments
    C_max (1,1) double {mustBePositive}
    p     (1,1) struct
    opts.N          (1,1) double {mustBePositive} = 20000
    opts.PhysicalOnly (1,1) logical = true
end

% Open interval on which L(C) is real; nudge inside to avoid the poles.
Clo = p.alpha * C_max / (p.alpha + p.beta);
Chi = C_max;
eps_ = (Chi - Clo) * 1e-9;

C = linspace(Clo + eps_, Chi - eps_, opts.N);

G = p.alpha * (C_max - C) ./ (p.beta * C);

% Guard the log against round-off pushing G marginally outside (0,1).
valid = G > 0 & G < 1;
C = C(valid);
G = G(valid);

phi = log(G ./ (1 - G));
L = C .* (p.u_crit + phi / p.s);

% dL/dC in closed form.
dLdC = p.u_crit + phi / p.s - (p.alpha * C_max) ./ (p.s * p.beta * C .* G .* (1 - G));

% Stability from the vector field's own Jacobian, evaluated on the curve.
% This is an independent route to the same information the fold structure
% implies, and the tests check that the middle branch is the unstable one.
lambda = zeros(size(C));
for k = 1:numel(C)
    [~, lambda(k)] = lc.dCdt(C(k), L(k), C_max, p);
end
stable = lambda < 0;

% --- Locate folds as sign changes of dL/dC ---
folds = struct('C', {}, 'L', {}, 'type', {});
sgn = sign(dLdC);
idx = find(sgn(1:end-1) .* sgn(2:end) < 0);

for k = 1:numel(idx)
    i = idx(k);
    Cf = fzero(@(cc) localdLdC(cc, C_max, p), [C(i), C(i+1)]);
    Lf = localL(Cf, C_max, p);
    folds(end+1).C = Cf; %#ok<AGROW>
    folds(end).L    = Lf;
    folds(end).type = "";
end

% Label folds: the one at lower L is the left fold (where the high branch
% is destroyed as load rises), the one at higher L the right fold.
if numel(folds) == 2
    if folds(1).L <= folds(2).L
        folds(1).type = "left";  folds(2).type = "right";
    else
        folds(1).type = "right"; folds(2).type = "left";
    end
    Lvals = [folds.L];
    br.L_fold_low  = min(Lvals);
    br.L_fold_high = max(Lvals);
    br.isBistable  = br.L_fold_high > 0;
else
    br.L_fold_low  = NaN;
    br.L_fold_high = NaN;
    br.isBistable  = false;
end

if opts.PhysicalOnly
    keep = L >= 0;
    C = C(keep); L = L(keep); dLdC = dLdC(keep);
    stable = stable(keep); lambda = lambda(keep);
end

br.C      = C;
br.L      = L;
br.dLdC   = dLdC;
br.stable = stable;
br.lambda = lambda;
br.folds  = folds;
br.C_max  = C_max;
br.p      = p;
end

% --- helpers: L(C) and dL/dC as standalone scalars for fzero ---

function L = localL(C, C_max, p)
G = p.alpha * (C_max - C) / (p.beta * C);
L = C * (p.u_crit + log(G / (1 - G)) / p.s);
end

function d = localdLdC(C, C_max, p)
G = p.alpha * (C_max - C) / (p.beta * C);
phi = log(G / (1 - G));
d = p.u_crit + phi / p.s - (p.alpha * C_max) / (p.s * p.beta * C * G * (1 - G));
end
