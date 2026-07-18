function [C0, info] = initialCapacity(p, isSeptic)
%INITIALCAPACITY Starting capacity, with sepsis applied where the data put it.
%
%   [C0, info] = vidd.initialCapacity(p, isSeptic) returns the capacity at
%   ICU admission.
%
%   SEPSIS IS AN OFFSET, NOT A RATE -- AND THAT IS THE MODEL'S BIGGEST
%   CORRECTION AGAINST ITS OWN SPEC
%   ------------------------------------------------------------------
%   The build plan puts sepsis in the dynamics, as d_disease*C: a persistent
%   activity-independent drain. Two studies say otherwise, and they are the
%   best human data available:
%
%     Demoule 2013 (PMID 23641946, n=85) -- sepsis was an independent
%       predictor of low twitch pressure ALREADY AT ADMISSION (coefficient
%       -3.74 cmH2O, P=0.002), and "Day 1 and Day 3 Ptr,stim were similar".
%       No measurable decline over two days. That is an initial condition,
%       not a slope.
%
%     Lecronier 2022 (PMID 35403916, n=92) -- septic patients' twitch
%       pressure ROSE 19% over a median 4 days while still ventilated,
%       against -7% in non-septic patients (p=0.005). Sepsis-associated
%       diaphragm dysfunction is reversible; the paper's title says so.
%
%   So under p.d_mode = 'offset' (the default) sepsis reduces C(0) and
%   d_disease is zero. Under 'rate' the spec's drain is restored, C(0) is
%   untouched, and the two can be compared -- see
%   experiments/vidd_exp3_catabolic.m.
%
%   THE OFFSET IS APPLIED FRACTIONALLY, AND IT HAS TO BE
%   ---------------------------------------------------
%   Demoule's -3.74 cmH2O is on the TWITCH PRESSURE scale, where a
%   non-septic patient sits at ~10 cmH2O. Subtracting 3.74 cmH2O from an MIP
%   of 80 would be a category error -- it would transfer a 36% effect as a
%   5% one. It is applied as the fractional reduction it represents, which
%   is the same declared proportional-loss assumption that lets any of this
%   literature reach the MIP scale at all.
%
%   See also vidd.simulateCapacity, vidd.dCdt

arguments
    p        (1,1) struct
    isSeptic (1,1) logical = false
end

info.isSeptic = isSeptic;
info.d_mode   = p.d_mode;
info.C_max0   = p.C_max0;

switch lower(string(p.d_mode))
    case "offset"
        if isSeptic
            C0 = p.C_max0 * (1 - p.sepsis_offset_fraction);
            info.appliedFraction = p.sepsis_offset_fraction;
        else
            C0 = p.C_max0;
            info.appliedFraction = 0;
        end
        info.d_disease_used = 0;

    case "rate"
        % The spec's version: sepsis is a drain, so the patient starts
        % healthy and declines. Retained for the comparison, not because the
        % data support it.
        C0 = p.C_max0;
        info.appliedFraction = 0;
        info.d_disease_used = p.d_disease;

    otherwise
        error('vidd:initialCapacity:unknownMode', ...
            'd_mode must be ''offset'' or ''rate'' (got "%s").', p.d_mode);
end

info.C0 = C0;
end
