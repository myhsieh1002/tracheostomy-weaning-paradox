function [K1, K2, fitInfo] = calibrateRohrer(flow, dP)
%CALIBRATEROHRER Least-squares fit of Rohrer coefficients to bench data.
%
%   [K1, K2] = wob.calibrateRohrer(flow, dP) fits
%
%       dP = K1*flow + K2*flow*|flow|
%
%   to measured (flow [L/s], dP [cmH2O]) pairs. The model is linear in the
%   coefficients, so this is an ordinary linear least-squares solve with no
%   starting guess and no iteration.
%
%   [K1, K2, fitInfo] = ... also returns a struct with the R-squared, RMSE,
%   residuals, and the number of points used.
%
%   This is the interface through which published bench data replaces the
%   placeholder coefficients in config/params_m1.json. Pass digitised
%   (flow, dP) points from a paper and write the returned K1/K2 back into
%   the config together with the citation.
%
%   Example:
%       % Bench points read off a published pressure-drop curve
%       flow = [0.25 0.5 0.75 1.0 1.25];      % L/s
%       dP   = [1.8 4.4 7.8 12.0 17.0];       % cmH2O
%       [K1, K2, info] = wob.calibrateRohrer(flow, dP);
%
%   See also wob.rohrerDrop

arguments
    flow (:,1) double {mustBeFinite}
    dP   (:,1) double {mustBeFinite}
end

if numel(flow) ~= numel(dP)
    error('wob:calibrateRohrer:sizeMismatch', ...
        'flow and dP must have the same number of elements (got %d and %d).', ...
        numel(flow), numel(dP));
end
if numel(flow) < 2
    error('wob:calibrateRohrer:tooFewPoints', ...
        'Need at least 2 points to fit 2 coefficients (got %d).', numel(flow));
end

% Design matrix: columns are the linear and quadratic Rohrer terms.
A = [flow, flow .* abs(flow)];
coeffs = A \ dP;

K1 = coeffs(1);
K2 = coeffs(2);

if nargout > 2
    predicted = A * coeffs;
    residuals = dP - predicted;
    ssRes = sum(residuals.^2);
    ssTot = sum((dP - mean(dP)).^2);

    fitInfo.R2        = 1 - ssRes / ssTot;
    fitInfo.RMSE      = sqrt(mean(residuals.^2));
    fitInfo.residuals = residuals;
    fitInfo.predicted = predicted;
    fitInfo.nPoints   = numel(flow);
end
end
