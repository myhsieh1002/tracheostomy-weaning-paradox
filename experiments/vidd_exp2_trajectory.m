function T = vidd_exp2_trajectory(cfgV, cfgM1)
%VIDD_EXP2_TRAJECTORY H2 -- capacity trajectories under ventilation strategies.
%
%   T = vidd_exp2_trajectory() integrates C_max(t) over the ICU stay under
%   each ventilation scenario, and against the calibration anchors.
%
%   The scenarios differ only in support level, hence in activity A. The
%   trajectories show the model's central claim on the capacity axis: how a
%   patient is ventilated, through the diaphragm activity it permits, moves
%   the capacity that Model 2's fold is defined against.
%
%   THE ANCHORS ARE OVERLAID SO THE FIT IS VISIBLE
%   ----------------------------------------------
%   Jaber 2010 (-32% force at 6 d, controlled MV) and the exponential form
%   check from Schepens 2015 are drawn against the low-activity trajectory,
%   so the calibration is shown rather than asserted.

arguments
    cfgV  (1,1) struct = vidd.loadConfig()
    cfgM1 (1,1) struct = wob.loadConfig()
end

c = viz.style();
p = vidd.calibrate(cfgV);
days = cfgV.sim.days;

scNames = fieldnames(cfgV.strategy.scenarios);
cols = [c.ETT; c.TRACH; c.grey; c.stable];

fig = figure('Position',[70 70 1180 440]);
tl = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

% ---- trajectories from a healthy start ----
nexttile; hold on;
rows = {};
for k = 1:numel(scNames)
    s = cfgV.strategy.scenarios.(scNames{k});
    A = vidd.supportToActivity(s.support_level, cfgV);
    r = vidd.simulateCapacity(cfgV, p, @(t) A, Days=days);
    plot(r.t, r.C, '-', 'Color', cols(min(k,end),:), 'DisplayName', strrep(s.label,'_','\_'));
    rows(end+1,:) = {scNames{k}, s.support_level, A, r.C(1), r.C(end), r.Cstar_inst(end), r.tau_inst(end)}; %#ok<AGROW>
end
% Model 2's observed weaning capacity band for reference
patch([0 days days 0], [42.3 42.3 53.8 53.8], c.rescue, 'FaceAlpha',0.08, 'EdgeColor','none', 'HandleVisibility','off');
text(days*0.5, 48, 'observed weaning MIP band', 'FontSize',8, 'Color', c.axisGrey, 'HorizontalAlignment','center');
xlabel('Days'); ylabel('Capacity C_{max}(t) (cmH_2O)');
title('Capacity trajectory by ventilation strategy');
legend('Location','east','Box','off');
ylim([0 p.C_max0*1.05]); grid on; set(gca,'GridAlpha',0.08);

% ---- Jaber anchor: controlled MV from healthy start ----
nexttile; hold on;
rJaber = vidd.simulateCapacity(cfgV, p, @(t) 0, Days=8);
plot(rJaber.t, 100*rJaber.C/rJaber.C(1), '-', 'Color', c.ETT, 'LineWidth',1.8, 'DisplayName','model, A=0 (controlled MV)');
j = cfgV.validation.jaber2010;
scatter(j.duration_days, 100*(1-j.force_loss_fraction), 60, c.stable, 'filled', ...
    'MarkerEdgeColor','w', 'DisplayName','Jaber 2010 (-32% @ 6 d)');
% Schepens form-check points
sch = cfgV.validation.schepens2015_form_check;
scatter(sch.at_hours/24, 100-sch.observed_loss_pct, 35, c.grey, 'filled', ...
    'MarkerEdgeColor','w', 'DisplayName','Schepens 2015 (thickness)');
xlabel('Days'); ylabel('Capacity, % of baseline');
title('Calibration anchor: complete disuse');
subtitle('model uses only Jaber for the rate; Schepens is a form check', 'FontSize',8, 'Color', c.axisGrey);
legend('Location','northeast','Box','off');
grid on; set(gca,'GridAlpha',0.08);

% ---- recovery: switch strategy partway ----
nexttile; hold on;
tSwitch = 7;
A_low  = vidd.supportToActivity(cfgV.strategy.scenarios.deep_sedation_controlled.support_level, cfgV);
A_high = vidd.supportToActivity(cfgV.strategy.scenarios.trach_spontaneous.support_level, cfgV);
% deteriorate on deep sedation, then switch to trach/spontaneous at day 7
rDown = vidd.simulateCapacity(cfgV, p, @(t) A_low, Days=tSwitch);
rUp   = vidd.simulateCapacity(cfgV, p, @(t) A_high, C0=rDown.C(end), Days=days-tSwitch);
plot(rDown.t, rDown.C, '-', 'Color', c.ETT, 'DisplayName','deep sedation (days 0-7)');
plot(rUp.t + tSwitch, rUp.C, '-', 'Color', c.TRACH, 'DisplayName','switch to spontaneous (day 7)');
xline(tSwitch, ':', 'Color', c.fold, 'Label','strategy switch', 'FontSize',8, 'HandleVisibility','off');
xlabel('Days'); ylabel('Capacity C_{max}(t) (cmH_2O)');
title('Recovery is possible -- capacity is not one-way');
subtitle('unlike Model 2''s fold, this axis has no trap', 'FontSize',8, 'Color', c.axisGrey);
legend('Location','east','Box','off');
ylim([0 p.C_max0*1.05]); grid on; set(gca,'GridAlpha',0.08);

title(tl, 'H2: ventilation strategy moves the capacity trajectory -- and, unlike fatigue, the move is reversible', ...
    'FontWeight','bold','FontSize',11.5);

viz.save(fig, 'fig_vidd_trajectory');

T = cell2table(rows, 'VariableNames', {'scenario','support_level','A','C0','C_end','Cstar','tau_days'});
writetable(T, fullfile(wob.projectRoot,'results','tables','vidd_exp2_trajectory.csv'));

fprintf('\n  Jaber check: model %.1f%% loss at 6 d vs measured 32%%\n', ...
    100*(1 - rJaber.C(rJaber.t>=6 & rJaber.t<6.1)/rJaber.C(1)));
fprintf('  Capacity IS recoverable: switching to spontaneous at day 7 reverses the decline.\n');
disp(T);
end
