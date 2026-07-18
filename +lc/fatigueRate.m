function [g, dg_du] = fatigueRate(u, p)
%FATIGUERATE Sigmoidal fatigue activation as a function of utilisation.
%
%   g = lc.fatigueRate(u, p) evaluates
%
%       g(u) = 1 / (1 + exp(-s*(u - u_crit)))
%
%   the fraction of the maximal fatigue rate engaged at utilisation u.
%   g rises steeply as u crosses u_crit, which is what makes the positive
%   feedback -- and hence the bistability -- possible.
%
%   [g, dg_du] = ... also returns the derivative, g' = s*g*(1-g), used to
%   build the analytic Jacobian.
%
%   NOTE ON u_crit: u must be the ratio of MEAN INSPIRATORY pressure to
%   maximal inspiratory pressure, for which u_crit ~ 0.4 (Roussos & Macklem
%   1977). It is NOT the tension-time index, whose 0.15 threshold includes a
%   duty-cycle factor. See lc.utilisation, which enforces this pairing.
%
%   Implementation note: the logistic is evaluated in a branch-free stable
%   form. The naive 1/(1+exp(-x)) overflows for large negative x, which
%   arises routinely here because the continuation drives u far below u_crit
%   with s as large as 50.
%
%   See also lc.dCdt, lc.utilisation

arguments
    u double
    p (1,1) struct
end

x = p.s * (u - p.u_crit);

% Stable logistic: use exp(x)/(1+exp(x)) where x < 0 so the exponential
% argument is never positive and cannot overflow.
g = zeros(size(x));
pos = x >= 0;
g(pos)  = 1 ./ (1 + exp(-x(pos)));
ex = exp(x(~pos));
g(~pos) = ex ./ (1 + ex);

if nargout > 1
    dg_du = p.s .* g .* (1 - g);
end
end
