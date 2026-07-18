function T = vidd_exp5_Agap(cfgV, cfgM1, cfgM2)
%VIDD_EXP5_AGAP The one unmeasured number, stress-tested.
%
%   T = vidd_exp5_Agap() sweeps the ETT->tracheostomy increase in diaphragm
%   activity (the "A-gap") and shows where the capacity channel overtakes
%   the load channel.
%
%   WHY THIS EXPERIMENT EXISTS
%   --------------------------
%   The programme's most striking number -- that a tracheostomy's effect on
%   muscle CAPACITY dwarfs its effect on breathing LOAD -- rests entirely on
%   one assumption: that a tracheostomy roughly doubles diaphragm activity
%   (A_ETT ~ 0.3 to A_TRACH ~ 0.6). A literature search found that number is
%   NOT anchored to data. Diaphragm activity has never been measured across
%   the ETT-vs-tracheostomy contrast (Link 2 below is unmeasured), and even
%   the sedation data that ARE strong do not reproduce a doubling.
%
%   THE EVIDENCE, HONESTLY
%   ----------------------
%   Link 1 (tracheostomy reduces sedation): SUPPORTED at RCT grade -- TracMan
%   5 vs 8 sedation-days, Meng meta WMD -6 days, holding under randomisation
%   so it is not merely a marker of the decision to wean.
%   Link 2 (less sedation raises diaphragm activity): PLAUSIBLE but never
%   measured in this comparison.
%
%   So the direction is solid and the magnitude is unknown. This experiment
%   therefore does not report a single answer; it reports the answer AS A
%   FUNCTION of the unmeasured gap, and marks where the conclusion flips.
%   That is the honest form of a result that hinges on one uncalibrated
%   number.

arguments
    cfgV  (1,1) struct = vidd.loadConfig()
    cfgM1 (1,1) struct = wob.loadConfig()
    cfgM2 (1,1) struct = lc.loadConfig()
end

c = viz.style();
p = vidd.calibrate(cfgV);
pM2 = lc.calibrate(cfgM2);
nf = lc.normalizedFolds(pM2);

% Load-axis effect is FIXED by Model 1 (measured, calibrated).
LE = lc.coupling(cfgM1, cfgM2, 'ETT_7_5').L_total;
LT = lc.coupling(cfgM1, cfgM2, 'TRACH_8_0').L_total;
dL = LE - LT;
loadEffect_Cmax = dL / nf.l_high;     % rescue-window width in C_max, cmH2O

% Capacity-axis effect DEPENDS on the unmeasured A-gap. Anchor the ETT end
% at a fixed activity and sweep the trach end upward.
A_ETT = vidd.supportToActivity(cfgV.strategy.scenarios.ett_partial_support.support_level, cfgV);
gaps = linspace(0, 0.45, 60);
capEffect = zeros(size(gaps));
for k = 1:numel(gaps)
    A_TRACH = min(A_ETT + gaps(k), 1);
    capEffect(k) = vidd.equilibriumCapacity(A_TRACH, p) - vidd.equilibriumCapacity(A_ETT, p);
end

% Where do the two channels cross?
crossIdx = find(capEffect >= loadEffect_Cmax, 1, 'first');
if isempty(crossIdx)
    gapCross = NaN;
else
    gapCross = interp1(capEffect, gaps, loadEffect_Cmax, 'linear');
end

% The scenario values currently assumed.
A_TRACH_assumed = vidd.supportToActivity(cfgV.strategy.scenarios.trach_spontaneous.support_level, cfgV);
gapAssumed = A_TRACH_assumed - A_ETT;

fig = figure('Position',[70 70 1180 430]);
tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

% ---- Panel A: the two channels vs the A-gap ----
nexttile; hold on;
plot(gaps, capEffect, '-', 'Color', c.TRACH, 'LineWidth',1.8, 'DisplayName','capacity channel (2b)');
plot(gaps, repmat(loadEffect_Cmax, size(gaps)), '-', 'Color', c.ETT, 'LineWidth',1.8, ...
    'DisplayName','load channel (Model 1, fixed)');
if ~isnan(gapCross)
    xline(gapCross, ':', 'Color', c.stable, 'Label',sprintf('channels cross at gap %.2f', gapCross), ...
        'FontSize',8, 'HandleVisibility','off');
end
xline(gapAssumed, '--', 'Color', c.grey, 'Label',sprintf('assumed gap %.2f', gapAssumed), ...
    'FontSize',8, 'HandleVisibility','off');
xlabel('A-gap = A_{TRACH} - A_{ETT}  (UNMEASURED)');
ylabel('effect on C_{max} axis (cmH_2O)');
title('Which channel dominates depends on the one number we cannot measure');
legend('Location','northwest','Box','off');
grid on; set(gca,'GridAlpha',0.08);

% ---- Panel B: ratio, with the evidence band ----
nexttile; hold on;
ratio = capEffect / loadEffect_Cmax;
plot(gaps, ratio, '-', 'Color', c.stable, 'LineWidth',1.8, 'HandleVisibility','off');
yline(1, ':', 'Color', c.axisGrey, 'Label','equal', 'FontSize',8, 'HandleVisibility','off');
% The sedation-proxy band: a "fraction not heavily sedated" proxy gives a
% SMALLER gap than the assumed doubling (ETT ~0.71 -> trach ~0.96 is ~0.25,
% but from a higher base; the defensible range is wide and low).
patch([0 0.15 0.15 0], [0 0 max(ratio)*1.05 max(ratio)*1.05], c.grey, ...
    'FaceAlpha',0.10, 'EdgeColor','none', 'DisplayName','sedation-proxy range (small gap)');
xline(gapAssumed, '--', 'Color', c.grey, 'HandleVisibility','off');
xlabel('A-gap = A_{TRACH} - A_{ETT}');
ylabel('capacity effect / load effect');
title('The 5x claim needs the large end of the gap');
subtitle('at small, sedation-plausible gaps the two channels are comparable', 'FontSize',8, 'Color', c.axisGrey);
legend('Location','northwest','Box','off');
grid on; set(gca,'GridAlpha',0.08);

title(tl, 'The capacity-vs-load verdict is a function of an unmeasured parameter, not a fixed result', ...
    'FontWeight','bold','FontSize',12);

viz.save(fig, 'fig_vidd_Agap');

% ---- Table ----
T = table(gaps', capEffect', repmat(loadEffect_Cmax,numel(gaps),1), ratio', ...
    'VariableNames', {'A_gap','capacity_effect_cmH2O','load_effect_cmH2O','ratio'});
writetable(T, fullfile(wob.projectRoot,'results','tables','vidd_exp5_Agap.csv'));

fprintf('\n  Load channel (Model 1, FIXED): dL=%.2f -> %.2f cmH2O of C_max\n', dL, loadEffect_Cmax);
fprintf('  Capacity channel (2b): DEPENDS on the A-gap, which is UNMEASURED.\n');
fprintf('    channels cross at A-gap = %.2f\n', gapCross);
fprintf('    at the assumed gap %.2f, capacity/load = %.1fx\n', gapAssumed, ...
    interp1(gaps, ratio, gapAssumed));
fprintf('    at a small (sedation-plausible) gap 0.10, capacity/load = %.1fx\n', ...
    interp1(gaps, ratio, 0.10));
fprintf('\n  CONCLUSION: the capacity channel overtakes the load channel only above an\n');
fprintf('  A-gap of ~%.2f. Whether it does is empirically open -- Link 2 (sedation->\n', gapCross);
fprintf('  activity) has never been measured across the ETT/trach contrast.\n');
end
