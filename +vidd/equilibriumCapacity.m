function [Cstar, tau] = equilibriumCapacity(A, p)
%EQUILIBRIUMCAPACITY Closed-form equilibrium capacity and time constant.
%
%   [Cstar, tau] = vidd.equilibriumCapacity(A, p) returns
%
%       syn   = k_syn*h(A)
%       loss  = k_deg*g(A) + d_disease
%       C*(A) = syn*C_max0 / (syn + loss)                  cmH2O
%       tau(A) = 1 / (syn + loss)                          days
%
%   Closed form: the system is linear in C (see vidd.dCdt), so setting
%   dC/dt = 0 solves directly. No root-finding, and the approach to C* is a
%   single exponential with time constant tau.
%
%   C* is a fraction syn/(syn + loss) of the ceiling C_max0 -- a tug of war
%   between rebuilding and losing, which is the physical reading of the
%   model and the reason the (C_max0 - C) factor belongs there (vidd.dCdt).
%
%   H1 -- THE INVERTED U -- IS A PROPERTY OF g ALONE
%   ------------------------------------------------
%   Set h0 = 1, so h is CONSTANT and the training effect is gone entirely.
%   Then C* = k_syn*C_max0/(k_syn + k_deg*g(A) + d_disease), which is STILL
%   an inverted U in A, peaking exactly where g bottoms out, at A*. The
%   shape is the degradation term's; it is not an artefact of tuning two
%   functions against each other, and tests/tModel2b.m asserts this.
%
%   With h0 < 1 the rising training effect pulls the peak to the RIGHT of
%   A*, and it can be pulled all the way to A = 1 if h's slope outweighs the
%   injury penalty -- at which point the model says the optimum is maximal
%   effort and contradicts the U-shape it was built to express. Where the
%   peak sits is therefore NOT a robust prediction; that a peak exists, and
%   that both extremes are worse, is.
%
%   H3 -- WHY ACTIVITY OPTIMISATION CANNOT RESCUE A CATABOLIC PATIENT
%   ----------------------------------------------------------------
%   d_disease enters `loss` and is INDEPENDENT of A. The best any strategy
%   can do is drive g(A) to its floor g0. As d_disease grows it dominates
%   the denominator, C* -> syn*C_max0/d_disease for every A, and the
%   ABSOLUTE gap between the best and worst achievable capacity collapses:
%   strategy stops buying cmH2O. That is H3, visible in the algebra rather
%   than only in a simulation -- vidd.strategyLeverage quantifies it, and
%   documents why the RELATIVE gap does not collapse the same way.
%
%   See also vidd.dCdt, vidd.strategyLeverage

arguments
    A double
    p (1,1) struct
end

h = vidd.synthesis(A, p);
g = vidd.degradation(A, p);

syn  = p.k_syn .* h;
loss = p.k_deg .* g + p.d_disease;

Cstar = syn .* p.C_max0 ./ (syn + loss);

if nargout > 1
    tau = 1 ./ (syn + loss);
end
end
