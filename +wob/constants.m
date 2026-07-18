function c = constants()
%CONSTANTS Physical constants and unit conversions used across Model 1.
%
%   All internal computation is done in cmH2O, L, s. Conversion to joules
%   happens only at the metric-reporting boundary.

c.CMH2O_L_TO_JOULE = 0.0981;   % 1 cmH2O*L = 0.0981 J
c.MMHG_TO_CMH2O    = 1.35951;  % 1 mmHg = 1.35951 cmH2O
c.ML_TO_L          = 1e-3;

% Gas properties at body temperature, saturated (BTPS) - used only by the
% physics-based ID-scaling path, never by the empirical Rohrer path.
c.AIR_DENSITY_KG_M3    = 1.075;    % ~37 degC, saturated
c.AIR_VISCOSITY_PA_S   = 1.89e-5;  % ~37 degC
end
