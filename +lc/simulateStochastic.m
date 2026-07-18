function out = simulateStochastic(t, L_t, C_max, p, opts)
%SIMULATESTOCHASTIC Euler-Maruyama integration of the noisy load-capacity system.
%
%   out = lc.simulateStochastic(t, L_t, C_max, p) integrates
%
%       dC = [alpha*(C_max - C) - beta*C*g(L/C)] dt + sigma*dW
%
%   along a prescribed load trajectory L_t, returning the capacity path.
%
%   Name-value options:
%       Sigma   noise amplitude as a FRACTION of C_max (default 0.02).
%               Scaling by C_max keeps the noise-to-signal ratio invariant
%               under the model's own scale invariance; a fixed absolute
%               sigma would mean something different at C_max = 20 than at
%               70, and the config value would silently change meaning
%               across the capacity grid.
%       C0      initial capacity (default: the high branch at L_t(1))
%       Seed    RNG seed (default 42)
%       Floor   lower clamp on C (default C_max*1e-3). u = L/C is undefined
%               at C = 0, and a noise excursion can otherwise push the state
%               non-positive between steps.
%
%   Fields of `out`: t, C, L, dt, sigma_abs, hitFloor.
%
%   See also lc.dCdt, lc.earlyWarning

arguments
    t     (1,:) double
    L_t   (1,:) double
    C_max (1,1) double {mustBePositive}
    p     (1,1) struct
    opts.Sigma (1,1) double {mustBeNonnegative} = 0.02
    opts.C0    (1,1) double = NaN
    opts.Seed  (1,1) double = 42
    opts.Floor (1,1) double = NaN
end

n = numel(t);
if numel(L_t) ~= n
    error('lc:simulateStochastic:sizeMismatch', ...
        't and L_t must have the same length (got %d and %d).', n, numel(L_t));
end

dt = diff(t);
if any(dt <= 0)
    error('lc:simulateStochastic:nonMonotonicTime', 't must be strictly increasing.');
end

sigma = opts.Sigma * C_max;
floorC = opts.Floor;
if isnan(floorC), floorC = C_max * 1e-3; end

% Start on the high branch unless told otherwise.
C0 = opts.C0;
if isnan(C0)
    fp = lc.fixedPoints(L_t(1), C_max, p);
    C0 = max(fp.C(fp.stable));
end

rng(opts.Seed);
dW = randn(1, n-1);

C = zeros(1, n);
C(1) = C0;
hitFloor = false;

for k = 1:n-1
    drift = lc.dCdt(C(k), L_t(k), C_max, p);
    Cnext = C(k) + drift*dt(k) + sigma*sqrt(dt(k))*dW(k);
    if Cnext < floorC
        Cnext = floorC;
        hitFloor = true;
    end
    C(k+1) = Cnext;
end

out.t         = t;
out.C         = C;
out.L         = L_t;
out.dt        = dt;
out.sigma_abs = sigma;
out.C_max     = C_max;
out.hitFloor  = hitFloor;
end
