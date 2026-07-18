function out = simulateCapacity(cfg, p, A_fun, opts)
%SIMULATECAPACITY Integrate the capacity trajectory C_max(t).
%
%   out = vidd.simulateCapacity(cfg, p, A_fun) integrates dC/dt over
%   cfg.sim.days with activity given by A_fun(t) (t in days).
%
%   out = vidd.simulateCapacity(..., C0=..., Days=...) overrides.
%
%   This is the function Model 2's spec section 2.5b names as the source of
%   C_max(t): its output feeds lc.rescueWindowTrajectory, which recomputes
%   the fold at each time point. This model does NOT re-derive the fold.
%
%   ode45 rather than a stiff solver: the system is linear in C with a
%   single time constant of order 1/(k_deg + d_disease) ~ 1-2 days, and
%   there is nothing fast to be stiff against. The stiffness in this
%   programme lives at the SEAM -- days here against hours in Model 2 -- and
%   is handled by the quasi-static coupling, not by the integrator.
%
%   Fields of `out`: t, C, A, Cstar_inst, tau_inst, p.
%
%   See also vidd.dCdt, vidd.supportToActivity, lc.rescueWindowTrajectory

arguments
    cfg   (1,1) struct
    p     (1,1) struct
    A_fun (1,1) function_handle
    opts.C0       (1,1) double  = NaN
    opts.Days     (1,1) double  = NaN
    opts.IsSeptic (1,1) logical = false
end

days = opts.Days;
if isnan(days), days = cfg.sim.days; end

% Sepsis enters through the INITIAL CONDITION by default, not through the
% dynamics -- see vidd.initialCapacity for why the data put it there.
C0 = opts.C0;
if isnan(C0)
    [C0, c0info] = vidd.initialCapacity(p, opts.IsSeptic);
else
    c0info = struct('C0', C0, 'isSeptic', opts.IsSeptic, 'note', 'C0 supplied explicitly');
end

t = 0:cfg.sim.dt_day:days;

odeOpts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
[tOut, C] = ode45(@(tt, CC) vidd.dCdt(CC, A_fun(tt), p), t, C0, odeOpts);

A = arrayfun(A_fun, tOut(:)');
[Cstar, tau] = vidd.equilibriumCapacity(A, p);

out.t          = tOut(:)';
out.C          = C(:)';
out.A          = A;
out.Cstar_inst = Cstar;    % the equilibrium the state is chasing at each t
out.tau_inst   = tau;      % days
out.p          = p;
out.C0         = C0;
out.c0info     = c0info;
out.isSeptic   = opts.IsSeptic;
end
