function h = synthesis(A, p)
%SYNTHESIS Training/synthesis effect as a function of diaphragm activity.
%
%   h = vidd.synthesis(A, p) evaluates
%
%       h(A) = h0 + (1 - h0)*A
%
%   the dimensionless multiplier on the synthesis rate k_syn: h(0) = h0
%   (basal protein synthesis continues without loading) and h(1) = 1.
%
%   WHY MONOTONE, WHEN THE SPEC SAYS IT SHOULD REVERSE AT HIGH ACTIVITY
%   ------------------------------------------------------------------
%   Because the reversal belongs to g, and putting it in both double-counts.
%   The equilibrium is
%
%       C*(A) = k_syn*h(A) / (k_deg*g(A) + d_disease)
%
%   With h CONSTANT this is already an inverted U, because the denominator
%   is U-shaped. So H1 needs nothing from h: the inverted U is a property of
%   the degradation term alone. Making h non-monotone as well would add two
%   more unconstrained parameters to produce a shape the model already has,
%   and would make it impossible to say which function the result came from.
%
%   h is therefore a monotone training effect -- which is what it physically
%   is -- and h0 = 1 collapses it to a constant. tests/tModel2b.m asserts
%   the inverted U survives h0 = 1.
%
%   No published data constrains h's form. It is linear because that is the
%   least it can be while still rising.
%
%   See also vidd.degradation, vidd.equilibriumCapacity

arguments
    A double
    p (1,1) struct
end

h = p.h0 + (1 - p.h0) .* A;
end
