function traj = rescueWindowTrajectory(t, C_max_t, L_device_range_t, p)
%RESCUEWINDOWTRAJECTORY Dynamic mode -- recompute the rescue window over time.
%
%   traj = lc.rescueWindowTrajectory(t, C_max_t, L_device_range_t, p) walks a
%   capacity trajectory and calls lc.rescueWindow at each time point,
%   producing the "does the window open or close over time" result
%   (Model 2's m2_exp6 / Model 2b's vidd_exp4 -- the same integration
%   figure).
%
%   INPUTS
%       t                 (1xN) time points (days)
%       C_max_t           (1xN) capacity ceiling at each time, from Model 2b
%       L_device_range_t  (Nx2) [L_TRACH, L_ETT] at each time, from Model 1
%                         (optionally via Model 1b). A 1x2 row is broadcast
%                         to every time point for the static-load case.
%       p                 dynamical parameters
%
%   STATUS: this is the INTEGRATION path. Model 2b does not exist yet, so
%   nothing in this repo currently supplies C_max_t. The function is
%   complete and tested against a synthetic decaying trajectory
%   (tests/tRescueApi.m) so that when 2b lands, the seam is already
%   specified and exercised: 2b calls lc.rescueWindow through here and does
%   NOT reimplement the fold solve.
%
%   Fields of `traj`: t, C_max, fold_left, fold_right, spans_fold,
%   window_lo, window_hi, windowOpen, firstClose.
%
%   See also lc.rescueWindow

arguments
    t                (1,:) double
    C_max_t          (1,:) double {mustBePositive}
    L_device_range_t (:,2) double {mustBeNonnegative}
    p                (1,1) struct
end

n = numel(t);
if numel(C_max_t) ~= n
    error('lc:rescueWindowTrajectory:sizeMismatch', ...
        't and C_max_t must have the same length (got %d and %d).', n, numel(C_max_t));
end

if size(L_device_range_t, 1) == 1
    L_device_range_t = repmat(L_device_range_t, n, 1);
elseif size(L_device_range_t, 1) ~= n
    error('lc:rescueWindowTrajectory:loadSizeMismatch', ...
        'L_device_range_t must have 1 or %d rows (got %d).', n, size(L_device_range_t,1));
end

fold_left  = nan(1,n); fold_right = nan(1,n);
spans_fold = false(1,n); windowOpen = false(1,n);
window_lo  = nan(1,n); window_hi  = nan(1,n);

for k = 1:n
    rw = lc.rescueWindow(C_max_t(k), L_device_range_t(k,:), p);
    fold_left(k)  = rw.fold_left;
    fold_right(k) = rw.fold_right;
    spans_fold(k) = rw.spans_fold;
    windowOpen(k) = ~isempty(rw.window);
    if ~isempty(rw.window)
        window_lo(k) = rw.window(1);
        window_hi(k) = rw.window(2);
    end
end

traj.t          = t;
traj.C_max      = C_max_t;
traj.fold_left  = fold_left;
traj.fold_right = fold_right;
traj.spans_fold = spans_fold;
traj.window_lo  = window_lo;
traj.window_hi  = window_hi;
traj.windowOpen = windowOpen;

% When does the window shut for good? The clinically load-bearing moment:
% after it, no device change can alter the outcome.
closedIdx = find(~windowOpen, 1, 'first');
if isempty(closedIdx)
    traj.firstClose = NaN;
else
    traj.firstClose = t(closedIdx);
end
end
