function pat = breathingPattern(pattern, opts)
%BREATHINGPATTERN Generate flow, volume and acceleration over one breath.
%
%   pat = wob.breathingPattern(pattern) builds a single respiratory cycle
%   from a struct with fields V_T_L, RR, IE_ratio, waveform and dt_s.
%
%   pat = wob.breathingPattern(pattern, ExpTimeConstant=tau) sets the
%   passive expiratory time constant in seconds (default 0.6 s).
%
%   Sign convention: inspiratory flow is POSITIVE, volume is measured above
%   FRC and decays towards zero during expiration.
%
%   Inspiratory waveforms (`pattern.waveform`), each scaled so that the
%   inspiratory flow integral is exactly V_T:
%       'sinusoidal'   A*sin(pi*t/Ti)      peak at mid-inspiration (default)
%       'half_sine'    alias of 'sinusoidal' -- for inspiratory flow the two
%                      names denote the same half-period sine. The build
%                      spec lists both; they are kept as aliases rather than
%                      invented as different shapes.
%       'constant'     A                   square / constant flow
%       'decelerating' A*(1 - t/Ti)        ramp down, peak at start
%       'accelerating' A*(t/Ti)            ramp up, peak at end
%
%   Expiration is modelled as PASSIVE exponential decay with time constant
%   tau, which is what happens on a T-piece or trach mask with no
%   ventilator support. Volume therefore approaches, but does not exactly
%   reach, zero within a finite expiratory time; the residual is returned
%   as pat.trappedVolume and is a genuine (if small) dynamic hyperinflation
%   signal rather than a numerical artefact.
%
%   Every Model 1 metric is computed over inspiration only, so the
%   expiratory limb does not enter WOB, PTP or f_device.
%
%   Fields of `pat`: t, flow, volume, accel, Ti, Te, Ttot, dutyCycle,
%   isInsp, trappedVolume, waveform, V_T, RR.
%
%   See also wob.simulateModeA, wob.simulateModeB

arguments
    pattern (1,1) struct
    opts.ExpTimeConstant (1,1) double {mustBePositive} = 0.6
end

V_T = pattern.V_T_L;
RR  = pattern.RR;
dt  = pattern.dt_s;

if V_T <= 0
    error('wob:breathingPattern:badTidalVolume', 'V_T_L must be positive (got %g).', V_T);
end
if RR <= 0
    error('wob:breathingPattern:badRate', 'RR must be positive (got %g).', RR);
end

Ttot = 60 / RR;
ie   = pattern.IE_ratio;
Ti   = Ttot * ie(1) / (ie(1) + ie(2));
Te   = Ttot - Ti;

% Build each limb on its own grid so that the inspiratory limb ENDS EXACTLY
% at Ti. A single 0:dt:Ttot grid truncates inspiration at the last sample
% below Ti -- for RR=18, I:E=1:2, dt=1 ms that is 1.111 s against a true Ti
% of 1.1111 s. Tidal volume still comes out right (the limb is normalised),
% but every flow-derived quantity is then off by ~1e-4 because the mean flow
% is computed over the wrong duration. The step size is adjusted by <1 part
% in 1e3 to make the grids commensurate, which is the cheaper error.
nIns = max(round(Ti/dt), 2) + 1;
nExp = max(round(Te/dt), 2) + 1;
tIns = linspace(0, Ti, nIns)';
tExpLocal = linspace(0, Te, nExp)';

t = [tIns; Ti + tExpLocal(2:end)];
isInsp = [true(nIns,1); false(nExp-1,1)];
tExp = tExpLocal(2:end);

waveform = lower(string(pattern.waveform));
switch waveform
    case {"sinusoidal", "half_sine", "half-sine"}
        flowIns = sin(pi * tIns / Ti);
    case {"constant", "square"}
        flowIns = ones(size(tIns));
    case {"decelerating", "decelerating_ramp"}
        flowIns = 1 - tIns / Ti;
    case {"accelerating", "accelerating_ramp"}
        flowIns = tIns / Ti;
    otherwise
        error('wob:breathingPattern:unknownWaveform', ...
            'Unknown waveform "%s". Valid: sinusoidal, half_sine, constant, decelerating, accelerating.', ...
            waveform);
end

% Scale the unit-shape inspiratory limb so that the DISCRETE trapezoidal
% integral equals V_T exactly. Normalising against the same quadrature rule
% that the metrics use keeps the analytic WOB checks exact rather than
% leaving them off by the quadrature error of the chosen dt.
flowIns = flowIns * (V_T / trapz(tIns, flowIns));

tau = opts.ExpTimeConstant;
flowExp = -(V_T / tau) * exp(-tExp / tau);

flow = [flowIns; flowExp];

% One cumulative integration over the whole cycle keeps volume continuous
% across the inspiratory-expiratory transition by construction.
volume = cumtrapz(t, flow);

% Flow acceleration, needed only when inertance is enabled. Expiratory
% onset is a step in flow, so differentiating across the transition would
% put a spurious spike on the last inspiratory sample. Each limb is
% therefore differentiated on its own.
accel = zeros(size(flow));
accel(isInsp)  = gradient(flowIns, tIns);
accel(~isInsp) = gradient(flowExp, tExp);

pat.t             = t;
pat.flow          = flow;
pat.volume        = volume;
pat.accel         = accel;
pat.Ti            = Ti;
pat.Te            = Te;
pat.Ttot          = Ttot;
pat.dutyCycle     = Ti / Ttot;
pat.isInsp        = isInsp;
pat.trappedVolume = volume(end);
pat.waveform      = waveform;
pat.V_T           = V_T;
pat.RR            = RR;
end
