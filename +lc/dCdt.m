function [f, dfdC] = dCdt(C, L, C_max, p)
%DCDT Load-capacity vector field and its analytic Jacobian.
%
%   f = lc.dCdt(C, L, C_max, p) evaluates
%
%       dC/dt = alpha*(C_max - C) - beta*C*g(L/C)
%
%   where C is sustainable inspiratory capacity (cmH2O), L is the
%   inspiratory load (cmH2O), and g is the sigmoidal fatigue activation.
%
%   [f, dfdC] = ... also returns df/dC, the 1-D Jacobian. A fixed point is
%   stable when dfdC < 0.
%
%   THE FEEDBACK
%   ------------
%   Recovery pulls C towards C_max. Fatigue removes capacity at a rate that
%   depends on utilisation u = L/C -- and u rises as C falls. Below u_crit
%   recovery dominates and the system sits near C_max; above it, falling C
%   raises u, which raises fatigue, which lowers C further. That runaway is
%   the positive feedback, and it is why sweeping L can produce one or three
%   fixed points and hence two saddle-node folds.
%
%   ANALYTIC JACOBIAN
%   -----------------
%   With u = L/C and du/dC = -L/C^2,
%
%       df/dC = -alpha - beta*g(u) + beta*s*g(u)*(1-g(u))*u
%
%   The last term is positive and is what can outweigh the two negative
%   terms to destabilise the middle branch. Supplying this in closed form
%   (rather than differencing) is what keeps the fold detection sharp: near
%   a fold df/dC passes through zero, and a finite-difference Jacobian there
%   is dominated by its own truncation error.
%
%   BEHAVIOUR AT C -> 0
%   -------------------
%   As C -> 0+, u -> inf and g -> 1, so f -> alpha*C_max > 0. C = 0 is
%   therefore NOT an attractor: the low branch is a low-but-positive
%   capacity, not collapse to zero. The "failure" state is a state of
%   sustained low capacity, which is the physiologically meaningful reading.
%
%   See also lc.fatigueRate, lc.fixedPoints

arguments
    C     double
    L     double
    C_max double
    p     (1,1) struct
end

if any(C <= 0, 'all')
    error('lc:dCdt:nonPositiveCapacity', ...
        'C must be strictly positive; u = L/C is undefined at C = 0.');
end

u = L ./ C;
[g, dg_du] = lc.fatigueRate(u, p);

f = p.alpha .* (C_max - C) - p.beta .* C .* g;

if nargout > 1
    % du/dC = -L/C^2, so  d/dC[ beta*C*g(u) ] = beta*g + beta*C*dg_du*(-L/C^2)
    %                                          = beta*g - beta*dg_du*u
    dfdC = -p.alpha - p.beta .* g + p.beta .* dg_du .* u;
end
end
