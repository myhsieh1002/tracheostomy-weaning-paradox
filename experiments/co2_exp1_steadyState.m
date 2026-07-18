function T = co2_exp1_steadyState(cfgCO2, cfgM1)
%CO2_EXP1_STEADYSTATE H1 -- steady-state PaCO2 vs ventilation, three arms.
%
%   T = co2_exp1_steadyState() sweeps minute ventilation and reports the
%   steady-state PaCO2 for the non-intubated, ETT and tracheostomy arms.
%
%   H1 IS A COUNTERFACTUAL AND THE FIGURE SAYS SO
%   ---------------------------------------------
%   "At matched V_E the tracheostomy gives a lower PaCO2" is true of the
%   alveolar gas equation and is not what happens to a patient. Chadda 2002
%   measured both states in the same subjects: removing 74 mL of dead space
%   left PaCO2 and respiratory rate UNCHANGED and raised V_T from 330 to
%   400 mL. PaCO2 is the regulated variable; ventilation is the response.
%
%   So this panel is the mechanism with the controller switched off. It is
%   worth showing because it isolates the dead-space effect, but the
%   clinical prediction is co2_exp2's, not this one's.

arguments
    cfgCO2 (1,1) struct = co2.loadConfig()
    cfgM1  (1,1) struct = wob.loadConfig()
end

c = viz.style();
phi = co2.alveolarDeadSpaceFraction(cfgCO2, cfgM1);

arms  = {cfgCO2.arms.native, cfgCO2.arms.ett, cfgCO2.arms.trach};
label = {'not intubated', 'ETT 7.5', 'trach 8.0'};
cols  = [c.grey; c.ETT; c.TRACH];

RR = cfgCO2.pattern.RR;
V_E = linspace(4, 16, 200);
V_T = V_E / RR;

rows = {};
fig = figure('Position',[70 70 1180 430]);
tl = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

% ---- Panel 1: PaCO2 vs V_E ----
nexttile; hold on;
P = nan(numel(arms), numel(V_E));
for a = 1:numel(arms)
    for k = 1:numel(V_E)
        try
            s = co2.steadyState(cfgCO2, cfgM1, arms{a}, V_T=V_T(k), RR=RR, Phi=phi);
            P(a,k) = s.PaCO2;
        catch ME
            if ~strcmp(ME.identifier, 'co2:steadyState:tidalBelowDeadSpace'), rethrow(ME); end
        end
    end
    plot(V_E, P(a,:), '-', 'Color', cols(a,:), 'DisplayName', label{a});
end
yline(cfgCO2.targets.PaCO2_target_mmHg, ':', 'Color', c.axisGrey, ...
    'Label','target PaCO_2', 'FontSize',8, 'HandleVisibility','off');
xlabel('Minute ventilation V_E (L/min)'); ylabel('Steady-state PaCO_2 (mmHg)');
title(sprintf('PaCO_2 at matched ventilation (RR = %g)', RR));
subtitle('a counterfactual: PaCO_2 is regulated, not V_E', 'FontSize',8, 'Color', c.axisGrey);
legend('Location','northeast','Box','off');
ylim([20 90]); grid on; set(gca,'GridAlpha',0.08);

% ---- Panel 2: PaCO2 vs VCO2 at fixed ventilation ----
nexttile; hold on;
VCO2s = linspace(150, 450, 200);
for a = 1:numel(arms)
    p = arrayfun(@(v) co2.steadyState(cfgCO2, cfgM1, arms{a}, VCO2=v, Phi=phi).PaCO2, VCO2s);
    plot(VCO2s, p, '-', 'Color', cols(a,:), 'DisplayName', label{a});
end
yline(cfgCO2.targets.PaCO2_target_mmHg, ':', 'Color', c.axisGrey, 'HandleVisibility','off');
xlabel('V̇CO_2 (mL/min)'); ylabel('Steady-state PaCO_2 (mmHg)');
title('PaCO_2 rises linearly with CO_2 production');
subtitle(sprintf('at the configured V_T = %.2f L, RR = %g', cfgCO2.pattern.V_T_L, RR), ...
    'FontSize',8, 'Color', c.axisGrey);
legend('Location','northwest','Box','off');
grid on; set(gca,'GridAlpha',0.08);

% ---- Panel 3: the gap between arms ----
nexttile; hold on;
gapBypass = P(1,:) - P(3,:);
gapDevice = P(2,:) - P(3,:);
plot(V_E, gapBypass, '-', 'Color', c.TRACH, 'DisplayName','not intubated \rightarrow trach');
plot(V_E, gapDevice, '-', 'Color', c.ETT,   'DisplayName','ETT \rightarrow trach');
xlabel('Minute ventilation V_E (L/min)'); ylabel('\DeltaPaCO_2 vs trach (mmHg)');
title('The dead-space gap, isolated');
subtitle('the device step is a fraction of the bypass step', 'FontSize',8, 'Color', c.axisGrey);
legend('Location','northeast','Box','off');
grid on; set(gca,'GridAlpha',0.08);

title(tl, 'H1: at matched ventilation, less dead space means a lower PaCO_2 — but ventilation is not what is held fixed in a patient', ...
    'FontWeight','bold','FontSize',11.5);

viz.save(fig, 'fig_co2_steady_state');

% ---- Table at the target ventilation ----
for a = 1:numel(arms)
    s = co2.steadyState(cfgCO2, cfgM1, arms{a}, Phi=phi);
    rows(end+1,:) = {label{a}, s.ds.total*1e3, s.ds.series*1e3, s.ds.alveolar*1e3, ...
        s.V_T, s.V_E, s.V_A, s.PaCO2}; %#ok<AGROW>
end
T = cell2table(rows, 'VariableNames', ...
    {'arm','Vd_total_mL','Vd_series_mL','Vd_alveolar_mL','V_T_L','V_E','V_A','PaCO2'});
writetable(T, fullfile(wob.projectRoot,'results','tables','co2_exp1_steady_state.csv'));

fprintf('\n  At the configured V_T = %.2f L, RR = %g, VCO2 = %g mL/min:\n', ...
    cfgCO2.pattern.V_T_L, cfgCO2.pattern.RR, cfgCO2.metabolism.VCO2_mL_min);
disp(T);
fprintf('  PaCO2 gap, not intubated -> trach: %.2f mmHg\n', T.PaCO2(1)-T.PaCO2(3));
fprintf('  PaCO2 gap, ETT -> trach:           %.2f mmHg  <- the device step\n', T.PaCO2(2)-T.PaCO2(3));
fprintf('\n  Chadda 2002 measured this contrast in vivo and found PaCO2 UNCHANGED:\n');
fprintf('  the patients raised V_T from 330 to 400 mL instead. This panel is the\n');
fprintf('  mechanism with the controller off; co2_exp2 is what the patient does.\n');
end
