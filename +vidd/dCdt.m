function [f, dfdC] = dCdt(C, A, p)
%DCDT VIDD capacity vector field and its analytic Jacobian.
%
%   f = vidd.dCdt(C, A, p) evaluates
%
%       dC/dt = k_syn*h(A)*(C_max0 - C) - k_deg*g(A)*C - d_disease*C
%
%   with C the sustainable inspiratory capacity (cmH2O, = Model 2's C_max)
%   and t in DAYS.
%
%   [f, dfdC] = ... also returns df/dC, which is constant in C.
%
%   THE (C_max0 - C) FACTOR IS NOT IN THE SPEC, AND HAS TO BE
%   --------------------------------------------------------
%   The build plan writes the synthesis term as k_syn*h(A) and labels k_syn
%   [1/day]. Those are inconsistent: for k_syn*h to be addable to dC/dt it
%   must carry cmH2O/day, not 1/day. Taken literally the equilibrium is
%   C* = k_syn*h/(k_deg*g + d_disease) = 0.4*0.6/0.7 = 0.34 cmH2O -- three
%   orders of magnitude below any physiological capacity.
%
%   Restoring the (C_max0 - C) factor fixes three things at once:
%     * k_syn's unit becomes 1/day, which is what the spec says it is;
%     * C_max0 enters the DYNAMICS instead of being an initial condition
%       that the spec lists among the parameters but never uses;
%     * the term reads as what it physically is -- the further below its
%       ceiling the muscle sits, the faster it rebuilds -- and mirrors
%       Model 2's own alpha*(C_max - C) recovery term.
%
%   The inverted U (H1) survives the change: see vidd.equilibriumCapacity.
%
%   THIS SYSTEM IS LINEAR IN C, AND THAT IS THE WHOLE DIFFERENCE FROM MODEL 2
%   ------------------------------------------------------------------------
%   Model 2's fatigue term is beta*C*g(L/C) -- the utilisation L/C puts C
%   inside the nonlinearity, which is what creates the positive feedback,
%   the multiple equilibria and the fold. Here every term is linear in C and
%   A is EXOGENOUS: no feedback from C, hence one equilibrium and no
%   bifurcation, at any parameter value.
%
%   So this model is monostable by construction. It does not tip; it drifts.
%   The tipping lives entirely in Model 2, and this model's job is to move
%   the parameter (C_max) that Model 2 tips with respect to. That division
%   is why the coupling is a clean slow-fast one rather than two dynamical
%   systems fighting each other, and tests/tModel2b.m asserts the
%   monostability rather than leaving it as an unstated hope.
%
%   THIS SYSTEM IS LINEAR IN C, AND THAT IS THE WHOLE DIFFERENCE FROM MODEL 2
%   ------------------------------------------------------------------------
%   Model 2's fatigue term is beta*C*g(L/C) -- the utilisation L/C puts C
%   inside the nonlinearity, which is what creates the positive feedback,
%   the multiple equilibria and the fold. Here the degradation is
%   k_deg*g(A)*C with A EXOGENOUS: no feedback from C, hence one equilibrium
%   and no bifurcation, at any parameter value.
%
%   So this model is monostable by construction. It does not tip; it drifts.
%   The tipping lives entirely in Model 2, and this model's job is to move
%   the parameter (C_max) that Model 2 tips with respect to. That division
%   is why the coupling is a clean slow-fast one rather than two dynamical
%   systems fighting each other, and tests/tModel2b.m asserts the
%   monostability rather than leaving it as an unstated hope.
%
%   See also vidd.degradation, vidd.synthesis, vidd.equilibriumCapacity

arguments
    C double
    A double
    p (1,1) struct
end

if any(A < 0 | A > 1, 'all')
    error('vidd:dCdt:activityOutOfRange', ...
        'A must lie in [0, 1]; got values in [%g, %g].', min(A(:)), max(A(:)));
end

h = vidd.synthesis(A, p);
g = vidd.degradation(A, p);

syn  = p.k_syn .* h;
loss = p.k_deg .* g + p.d_disease;

f = syn .* (p.C_max0 - C) - loss .* C;

if nargout > 1
    dfdC = -(syn + loss);
end
end
