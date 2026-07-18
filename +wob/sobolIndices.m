function S = sobolIndices(fun, ranges, opts)
%SOBOLINDICES First-order and total-effect Sobol indices (Saltelli estimators).
%
%   S = wob.sobolIndices(fun, ranges) estimates variance-based sensitivity
%   indices for a scalar model.
%
%   fun     function handle taking a 1xk row of parameter values, returning
%           a scalar
%   ranges  kx2 matrix of [lower, upper] bounds
%
%   Name-value options:
%       N          base sample size (default 1024). Total cost is N*(k+2).
%       Seed       RNG seed for reproducibility (default 42)
%       Names      1xk string array of parameter names
%       Bootstrap  bootstrap replicates for CIs (default 200; 0 to skip)
%
%   METHOD
%   ------
%   Two independent Sobol point sets A and B (scrambled, from `sobolset`),
%   plus k hybrid matrices AB_i = A with column i taken from B. Then
%   (Saltelli et al. 2010):
%
%       S1_i = mean( Y_B .* (Y_ABi - Y_A) ) / Var(Y)
%       ST_i = mean( (Y_A - Y_ABi).^2 ) / (2*Var(Y))
%
%   S1 is the variance explained by parameter i alone; ST additionally
%   counts every interaction it participates in. ST_i ~ S1_i means the
%   parameter acts independently; ST_i >> S1_i means it matters mainly
%   through interactions. Sum(S1) < 1 is the signature of interaction.
%
%   These estimators are used rather than the older Sobol/Jansen forms
%   because they are markedly less biased at small N, and N is the binding
%   cost here.
%
%   Fields of `S`: S1, ST, S1_ci, ST_ci, names, Y, varY, N, cost.
%
%   See also sobolset, wob.simulateEffort

arguments
    fun    (1,1) function_handle
    ranges (:,2) double
    opts.N         (1,1) double {mustBePositive} = 1024
    opts.Seed      (1,1) double = 42
    opts.Names     string = string.empty
    opts.Bootstrap (1,1) double {mustBeNonnegative} = 200
end

k = size(ranges, 1);
N = opts.N;

if any(ranges(:,2) <= ranges(:,1))
    bad = find(ranges(:,2) <= ranges(:,1), 1);
    error('wob:sobolIndices:badRange', ...
        'Range %d has upper <= lower ([%g, %g]).', bad, ranges(bad,1), ranges(bad,2));
end

if isempty(opts.Names)
    opts.Names = "p" + string(1:k);
end

% Two independent sample sets. One 2k-dimensional Sobol sequence split in
% half gives A and B the low-discrepancy property jointly, which is what the
% Saltelli scheme assumes.
sob = sobolset(2*k, 'Skip', 1, 'Leap', 0);
sob = scramble(sob, 'MatousekAffineOwen');
sob.Skip = 1;
rng(opts.Seed);
X = net(sob, N);

A = scaleTo(X(:, 1:k),      ranges);
B = scaleTo(X(:, k+1:2*k),  ranges);

Y_A = evalAll(fun, A);
Y_B = evalAll(fun, B);

Y_AB = zeros(N, k);
for i = 1:k
    AB = A;
    AB(:, i) = B(:, i);
    Y_AB(:, i) = evalAll(fun, AB);
end

[S1, ST, varY] = estimate(Y_A, Y_B, Y_AB);

S.S1    = S1;
S.ST    = ST;
S.varY  = varY;
S.names = opts.Names;
S.Y     = [Y_A; Y_B];
S.N     = N;
S.cost  = N * (k + 2);

% Bootstrap CIs over the sample index, which is the only randomness that
% matters once the design is fixed.
if opts.Bootstrap > 0
    nb = opts.Bootstrap;
    S1b = zeros(nb, k); STb = zeros(nb, k);
    for b = 1:nb
        idx = randi(N, N, 1);
        [S1b(b,:), STb(b,:)] = estimate(Y_A(idx), Y_B(idx), Y_AB(idx,:));
    end
    S.S1_ci = prctile(S1b, [2.5 97.5], 1);
    S.ST_ci = prctile(STb, [2.5 97.5], 1);
else
    S.S1_ci = nan(2, k);
    S.ST_ci = nan(2, k);
end
end

function [S1, ST, varY] = estimate(Y_A, Y_B, Y_AB)
k = size(Y_AB, 2);
varY = var([Y_A; Y_B], 1);

S1 = zeros(1, k); ST = zeros(1, k);
if varY == 0
    return;   % a constant model has no variance to apportion
end
for i = 1:k
    S1(i) = mean(Y_B .* (Y_AB(:,i) - Y_A)) / varY;
    ST(i) = mean((Y_A - Y_AB(:,i)).^2) / (2 * varY);
end
end

function P = scaleTo(U, ranges)
P = ranges(:,1)' + U .* (ranges(:,2)' - ranges(:,1)');
end

function Y = evalAll(fun, P)
n = size(P, 1);
Y = zeros(n, 1);
for j = 1:n
    Y(j) = fun(P(j, :));
end
end
