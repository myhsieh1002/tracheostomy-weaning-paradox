function S = vidd_exp6_calibrationSensitivity(cfgV, cfgM1, cfgM2)
%VIDD_EXP6_CALIBRATIONSENSITIVITY How much do the two anchor readings move Model 2b?
%
%   S = vidd_exp6_calibrationSensitivity() re-solves the Model 2b calibration
%   under the two readings of the literature that a 2026-11-18 re-check of
%   the primary sources left open, and reports what moves.
%
%   (1) JABER'S INTERVAL -> k_deg
%       Jaber 2011 reports -32% twitch pressure "after 5-6 days" of MV, with
%       the baseline taken a mean 12.5 h after intubation. The default reads
%       this as 5.5 d (k_deg = 0.0701/day); 5 d (0.0771) and 6 d (0.0643, the
%       v1.0.0 value) bracket it. k_deg feeds the Lecronier two-point solve,
%       so k_syn, the implied septic offset, and the capacity channel's
%       A-gap crossover all move with it.
%
%   (2) DEMOULE'S SEPTIC OFFSET -> sepsis_offset_fraction
%       Demoule 2013 reports only a regression coefficient (-3.74 cmH2O,
%       multivariable) and whole-population Day-1 Ptr,stim (median 8.2, mean
%       10.0); the non-septic reference has to be reconstructed, and the
%       fraction is 0.30-0.36 depending on how. The offset changes only C(0)
%       for septic scenarios and the comparison against Lecronier; it does
%       not enter the equilibria, so it does not move the A-gap crossover.
%
%   Nothing here changes the default calibration: it is reported alongside.

arguments
    cfgV  (1,1) struct = vidd.loadConfig()
    cfgM1 (1,1) struct = wob.loadConfig()
    cfgM2 (1,1) struct = lc.loadConfig()
end

j     = cfgV.validation.jaber2010;
days  = unique([j.duration_days_sensitivity(:)', j.duration_days]);
offs  = cfgV.disease.sepsis_offset_sensitivity(:)';

% Load channel (fixed by Model 1), as in vidd_exp5.
pM2 = lc.calibrate(cfgM2);
nf  = lc.normalizedFolds(pM2);
dL  = lc.coupling(cfgM1, cfgM2, 'ETT_7_5').L_total - lc.coupling(cfgM1, cfgM2, 'TRACH_8_0').L_total;
loadEffect = dL / nf.l_high;

A_ETT = vidd.supportToActivity(cfgV.strategy.scenarios.ett_partial_support.support_level, cfgV);
gaps  = linspace(0, 0.45, 400);

rows = zeros(numel(days), 8);
for i = 1:numel(days)
    c = cfgV;
    c.capacity.k_deg = -log(1 - j.force_loss_fraction) / days(i);
    c.validation.jaber2010.duration_days = days(i);
    [p, info] = vidd.calibrate(c);
    cap = arrayfun(@(g) vidd.equilibriumCapacity(min(A_ETT + g, 1), p), gaps) ...
          - vidd.equilibriumCapacity(A_ETT, p);
    ratio = cap / loadEffect;
    if any(cap >= loadEffect)
        gapCross = interp1(cap, gaps, loadEffect, 'linear');
    else
        gapCross = NaN;
    end
    rows(i,:) = [days(i), p.k_deg, p.k_syn, info.lecronier_septic_offset, ...
                 info.jaber_predicted_loss_at_A0, gapCross, ...
                 interp1(gaps, ratio, 0.10), interp1(gaps, ratio, 0.30)];
end
S.jaber = array2table(rows, 'VariableNames', {'interval_days','k_deg','k_syn', ...
    'lecronier_septic_offset','jaber_loss_reproduced','Agap_crossover', ...
    'ratio_at_gap_0_10','ratio_at_gap_0_30'});

% Offset sensitivity: against the DEFAULT-calibration Lecronier offset.
[~, info0] = vidd.calibrate(cfgV);
S.offset = table(offs', repmat(info0.lecronier_septic_offset, numel(offs), 1), ...
    abs(offs' - info0.lecronier_septic_offset), ...
    'VariableNames', {'demoule_offset','lecronier_septic_offset','gap'});
S.loadEffect = loadEffect;

outDir = fullfile(wob.projectRoot, 'results', 'tables');
writetable(S.jaber,  fullfile(outDir, 'vidd_exp6_jaber_interval.csv'));
writetable(S.offset, fullfile(outDir, 'vidd_exp6_sepsis_offset.csv'));

fprintf('\n  Jaber interval sensitivity (k_deg = -ln(0.68)/days):\n');
disp(S.jaber);
fprintf('  Demoule offset sensitivity (Lecronier-implied offset %.3f at default):\n', ...
    info0.lecronier_septic_offset);
disp(S.offset);
end
