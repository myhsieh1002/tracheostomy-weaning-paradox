function r = rescueOutcome(L_ETT, L_TRACH, C_max, p)
%RESCUEOUTCOME Does switching ETT -> tracheostomy flip the weaning outcome?
%
%   r = lc.rescueOutcome(L_ETT, L_TRACH, C_max, p) classifies one patient:
%   given the load under each device and the capacity ceiling, does the
%   tracheostomy move the operating point across the right fold?
%
%   THE CLOSED FORM
%   ---------------
%   By the scale invariance, a sustainable branch exists iff the normalised
%   load sits below the right fold:
%
%       L / C_max  <  l_high        <=>        L  <  l_high * C_max
%
%   So a tracheostomy rescues the patient exactly when the fold falls
%   BETWEEN the two loads:
%
%       L_TRACH  <  l_high * C_max  <  L_ETT
%
%   Rearranged, rescue occurs iff
%
%       C_max  in  ( L_TRACH / l_high ,  L_ETT / l_high )
%
%   -- a band in C_max whose WIDTH is (L_ETT - L_TRACH) / l_high, i.e.
%   directly proportional to the device load difference. This is the
%   rescue window, and it is analytic: no simulation is needed to find it.
%
%   THE PARADOX, STATED IN ONE LINE
%   -------------------------------
%   Below the band, both devices leave the patient above the fold and the
%   tracheostomy changes nothing. Above it, both sit below the fold and the
%   patient would have weaned regardless. Only inside the band does the
%   device decide the outcome -- and the band is only as wide as the load
%   the device actually removes.
%
%   Outcome codes in r.outcome:
%       "rescued"        ETT fails, trach weans      <- the window
%       "both_wean"      neither device is limiting
%       "both_fail"      capacity too low for either
%       "paradoxical"    trach fails where ETT weans (should be impossible
%                        for L_TRACH < L_ETT; flagged, not silently ignored)
%
%   See also lc.normalizedFolds, lc.coupling

arguments
    L_ETT   (1,1) double
    L_TRACH (1,1) double
    C_max   (1,1) double {mustBePositive}
    p       (1,1) struct
end

nf = lc.normalizedFolds(p);

if ~nf.isBistable
    error('lc:rescueOutcome:noBistability', ...
        ['The parameter set is monostable, so no fold exists and "rescue" is ' ...
         'undefined. Check beta/alpha and s.']);
end

L_fold = nf.l_high * C_max;

ettWeans   = L_ETT   < L_fold;
trachWeans = L_TRACH < L_fold;

if ~ettWeans && trachWeans
    outcome = "rescued";
elseif ettWeans && trachWeans
    outcome = "both_wean";
elseif ~ettWeans && ~trachWeans
    outcome = "both_fail";
else
    outcome = "paradoxical";
end

r.outcome      = outcome;
r.isRescued    = outcome == "rescued";
r.L_fold       = L_fold;
r.l_high       = nf.l_high;
r.l_ETT        = L_ETT / C_max;
r.l_TRACH      = L_TRACH / C_max;
r.ettWeans     = ettWeans;
r.trachWeans   = trachWeans;
r.C_max        = C_max;
r.dL_device    = L_ETT - L_TRACH;

% The analytic window this patient's loads imply.
r.window_C_max = [L_TRACH / nf.l_high, L_ETT / nf.l_high];
r.window_width = r.dL_device / nf.l_high;
end
