function A = supportToActivity(supportLevel, cfg)
%SUPPORTTOACTIVITY Map ventilator support onto diaphragm activity.
%
%   A = vidd.supportToActivity(supportLevel, cfg) evaluates
%
%       A = A_max * (1 - support_level)
%
%   with support_level 1 = fully controlled ventilation (the diaphragm does
%   nothing) and 0 = unsupported (it does everything).
%
%   THIS MAPPING IS THE MODEL'S WEAKEST LINK, AND IT IS LOAD-BEARING
%   ---------------------------------------------------------------
%   "Support level" is not a measured quantity. What the VIDD literature
%   actually measures is the diaphragm thickening fraction (TFdi) -- and the
%   empirical optimum, Goligher 2018's TFdi of 15-30%, is expressed in that
%   unit, not in this one. Treating support level as a linear proxy for
%   activity, and A* as commensurate with a TFdi band, is an assumption with
%   no calibration behind it.
%
%   Everything downstream that depends on WHERE a scenario sits relative to
%   A* inherits that assumption. The qualitative claims -- that there is an
%   optimum, that both extremes are worse, that catabolism erodes the
%   leverage of choosing well -- do not depend on the mapping being exactly
%   right. The claim that a particular support level IS optimal does, and
%   should not be made.
%
%   See also vidd.degradation

arguments
    supportLevel double
    cfg (1,1) struct
end

if any(supportLevel < 0 | supportLevel > 1, 'all')
    error('vidd:supportToActivity:outOfRange', ...
        'support_level must lie in [0, 1]; got values in [%g, %g].', ...
        min(supportLevel(:)), max(supportLevel(:)));
end

A = cfg.strategy.A_max * (1 - supportLevel);
end
