function out = physicsScaling(ID_mm, length_cm, flow)
%PHYSICSSCALING Physics-based pressure drop for a smooth circular tube.
%
%   out = wob.physicsScaling(ID_mm, length_cm, flow) returns a struct with
%   the laminar (Poiseuille) and turbulent (Blasius) pressure-drop
%   components in cmH2O for flow in L/s.
%
%   This is the AUXILIARY path of the model. It exists to show how
%   resistance scales with internal diameter -- the check behind H1 -- and
%   is deliberately NOT used to set the working K1/K2, which come from
%   published bench data via wob.calibrateRohrer. A real tube has
%   connectors, curvature and a bevel; this idealisation would understate
%   its resistance.
%
%   Scaling laws returned:
%       laminar   dP = 128*mu*L*Q / (pi*d^4)            ~ 1/d^4
%       turbulent dP = 0.241*rho^0.75*mu^0.25*L*Q^1.75 / d^4.75
%                                                        ~ 1/d^4.75
%   The turbulent expression is the Blasius friction factor
%   (f = 0.316*Re^-0.25) substituted into the Darcy-Weisbach equation.
%
%   Fields of `out`: dP_laminar, dP_turbulent, dP_total, Re, regime.

arguments
    ID_mm     (1,1) double {mustBePositive}
    length_cm (1,1) double {mustBePositive}
    flow      double
end

c = wob.constants();
rho = c.AIR_DENSITY_KG_M3;
mu  = c.AIR_VISCOSITY_PA_S;

d = ID_mm * 1e-3;         % m
L = length_cm * 1e-2;     % m
Q = abs(flow) * 1e-3;     % L/s -> m^3/s

area = pi * d^2 / 4;
velocity = Q ./ area;
Re = rho .* velocity .* d ./ mu;

% Poiseuille (Pa)
dP_lam_Pa = 128 * mu * L .* Q ./ (pi * d^4);

% Darcy-Weisbach with Blasius friction factor (Pa). Guard Re -> 0, where
% the Blasius correlation is undefined and the turbulent term is physically
% absent anyway.
fDarcy = zeros(size(Re));
nonzero = Re > 0;
fDarcy(nonzero) = 0.316 .* Re(nonzero).^(-0.25);
dP_turb_Pa = fDarcy .* (L / d) .* (rho / 2) .* velocity.^2;

PA_TO_CMH2O = 1 / 98.0665;

out.dP_laminar   = dP_lam_Pa * PA_TO_CMH2O .* sign(flow);
out.dP_turbulent = dP_turb_Pa * PA_TO_CMH2O .* sign(flow);
out.dP_total     = out.dP_laminar + out.dP_turbulent;
out.Re           = Re;
out.regime       = repmat("laminar", size(Re));
out.regime(Re > 2300) = "turbulent";
end
