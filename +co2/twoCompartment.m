function out = twoCompartment(cfgCO2, tSpan, V_A_fun, opts)
%TWOCOMPARTMENT Lung + body CO2 kinetics (the transient path).
%
%   out = co2.twoCompartment(cfgCO2, tSpan, V_A_fun) integrates
%
%       S_lung * dPA/dt = Q*slope*(Pv - PA)*1000  -  V_A*PA/K
%       S_body * dPv/dt = VCO2                    -  Q*slope*(Pv - PA)*1000
%
%   where V_A_fun(t) returns alveolar ventilation (L/min) at time t (min).
%
%   THE UNIT CHAIN, WHICH IS WHERE THIS MODEL IS EASIEST TO GET WRONG
%   -----------------------------------------------------------------
%   VCO2 is quoted STPD; alveolar gas is BTPS. Rather than carry an explicit
%   conversion, both elimination terms are written through the same constant
%   K = 0.863 that defines the alveolar gas equation, which already contains
%   it:
%
%       PA = K*VCO2/V_A     =>     VCO2_eliminated = V_A*PA/K
%
%   so the steady state of this ODE reproduces the algebraic relation by
%   construction rather than by coincidence. The lung capacitance follows the
%   same route: S_lung = V_L * (STPD/BTPS) / (PB - PH2O), in mL CO2(STPD) per
%   mmHg.
%
%   Transport is Q*slope*(Pv - PA), with slope the linear CO2 dissociation
%   coefficient. Only the SLOPE appears -- the intercept cancels from the
%   venous-arterial difference -- which is why the crude linear dissociation
%   is tolerable here.
%
%   STIFFNESS
%   ---------
%   The lung time constant is ~35 s and the body's ~35 min: a ratio of ~60.
%   ode15s is used rather than ode45 for that reason; ode45 would take its
%   step from the lung and crawl through the body's timescale.
%
%   Name-value: VCO2, PA0, Pv0, and a Jacobian toggle.
%
%   Fields of `out`: t, PA, Pv, V_A, tau_lung, tau_body, S_lung, S_body.
%
%   See also co2.steadyState

arguments
    cfgCO2  (1,1) struct
    tSpan   (1,:) double
    V_A_fun (1,1) function_handle
    opts.VCO2 (1,1) double = NaN
    opts.PA0  (1,1) double = NaN
    opts.Pv0  (1,1) double = NaN
end

VCO2 = pick(opts.VCO2, cfgCO2.metabolism.VCO2_mL_min);
K    = cfgCO2.constants.alveolar_K;
PB   = cfgCO2.constants.PB_mmHg;
PH2O = cfgCO2.constants.PH2O_mmHg;

Q      = cfgCO2.stores.Q_L_min;
slope  = cfgCO2.stores.diss_slope;
V_L    = cfgCO2.stores.V_L_lung_mL;
S_body = cfgCO2.stores.S_body_mL_per_kg_per_mmHg * cfgCO2.stores.body_weight_kg;

% Lung capacitance in mL CO2 (STPD) per mmHg alveolar PCO2.
STPD_over_BTPS = (273/310) * ((PB - PH2O)/PB);
S_lung = V_L * STPD_over_BTPS / (PB - PH2O);

% Transport conductance: mL CO2/min per mmHg of venous-arterial difference.
Gt = Q * slope * 1000;

% Start from the steady state of the initial ventilation unless told otherwise.
VA0 = V_A_fun(tSpan(1));
PA0 = pick(opts.PA0, K * VCO2 / VA0);
Pv0 = pick(opts.Pv0, PA0 + VCO2 / Gt);

odeFun = @(t, y) [ (Gt*(y(2) - y(1)) - V_A_fun(t)*y(1)/K) / S_lung ;
                   (VCO2 - Gt*(y(2) - y(1)))              / S_body ];

% The system is linear in (PA, Pv) for a fixed V_A, so the Jacobian is
% exact and cheap; giving it to the stiff solver removes the finite
% differencing it would otherwise do at every step.
jac = @(t, y) [ (-Gt - V_A_fun(t)/K)/S_lung,  Gt/S_lung ;
                 Gt/S_body,                  -Gt/S_body ];

odeOpts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10, 'Jacobian', jac);
[t, y] = ode15s(odeFun, tSpan, [PA0; Pv0], odeOpts);

out.t        = t(:)';
out.PA       = y(:,1)';
out.Pv       = y(:,2)';
out.V_A      = arrayfun(V_A_fun, out.t);
out.S_lung   = S_lung;
out.S_body   = S_body;
out.Gt       = Gt;
out.VCO2     = VCO2;
out.tau_lung = S_lung / (VA0/K);      % min
out.tau_body = S_body / Gt;           % min
end

function v = pick(given, fallback)
if isnan(given), v = fallback; else, v = given; end
end
