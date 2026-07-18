function T = vidd_exp4_dynamicRescue(cfgV, cfgM1, cfgM2)
%VIDD_EXP4_DYNAMICRESCUE H4 -- the dynamic rescue window. PROGRAM CAPSTONE.
%
%   T = vidd_exp4_dynamicRescue() walks a capacity trajectory C_max(t) from
%   Model 2b through Model 2's fold, so the rescue window opens and closes
%   over time. This is Model 2's m2_exp6 and this model's vidd_exp4 -- the
%   same figure, and the point where all three models meet.
%
%   THE SEAM
%   --------
%   Model 2b supplies C_max(t) (vidd.simulateCapacity). Model 1 supplies the
%   loads L_ETT, L_TRACH (lc.coupling). Model 2's named API
%   (lc.rescueWindow, via lc.rescueWindowTrajectory) recomputes the fold at
%   each time point. Nothing here re-derives the fold -- Model 2's spec
%   section 9 forbids it, and vidd.couplingToModel2 enforces it.
%
%   THE DEVICE MOVES BOTH AXES, WHICH IS THE WHOLE PROGRAM IN ONE FIGURE
%   -------------------------------------------------------------------
%   A tracheostomy lowers the LOAD (Model 1) and, by enabling lighter
%   sedation, raises the CAPACITY trajectory (Model 2b). So the trach arm
%   both sits at a lower operating point AND climbs a rising C_max(t), while
%   the sedated-ETT arm sits higher on load and sinks. Whether the operating
%   point ends up below the fold -- weanable -- is the product of both.
%
%   ⚠ THE RESULT HINGES ON AN UNCALIBRATED ASSUMPTION
%   -------------------------------------------------
%   The capacity axis is driven entirely by the activity A each scenario
%   implies, and the support_level -> A mapping is NOT anchored to data (see
%   vidd.supportToActivity). The relative SIZE of the two axes -- currently
%   the capacity effect dominates the load effect several-fold -- is a
%   consequence of that mapping, not a measured result. The figure states
%   this on its face. What is robust: that the device acts on both axes and
%   that the capacity axis, being reversible, is not diluted by disease the
%   way the load axis is diluted by mechanics.

arguments
    cfgV  (1,1) struct = vidd.loadConfig()
    cfgM1 (1,1) struct = wob.loadConfig()
    cfgM2 (1,1) struct = lc.loadConfig()
end

c = viz.style();

% The contrast is only visible in a patient whose load sits NEAR the fold --
% otherwise both arms are weanable at every capacity and nothing moves. A
% low-compliance, high-metabolic patient (the population Model 2's rescue
% window actually reaches) puts the operating point in the interesting
% region. This is the same lesson as m2_exp3: a healthy-lung patient makes
% the device look irrelevant because the device was never the binding
% constraint there.
sick = struct('C_rs', 0.03, 'R_aw_native', 10);
cfgM1sick = cfgM1;
cfgM1sick.patient.target_VA_L_min = cfgM1.disease_grid.target_VA_L_min(end);  % high metabolic

% Two whole-patient scenarios: a tracheostomy patient (low load, spontaneous)
% vs an ETT patient kept sedated (higher load, low activity).
cpTrach = vidd.couplingToModel2(cfgV, cfgM1sick, cfgM2, 'trach_spontaneous',   PatientOverrides=sick);
cpETT   = vidd.couplingToModel2(cfgV, cfgM1sick, cfgM2, 'deep_sedation_controlled', PatientOverrides=sick);

pM2 = cpTrach.pM2;
nf  = lc.normalizedFolds(pM2);

fig = figure('Position',[50 50 1260 460]);
tl = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

% ---- Panel A: capacity trajectory vs the moving fold ----
nexttile; hold on;
% The fold, expressed as the C_max at which each device's load sits exactly
% on it: C_max_fold = L / l_high. Above this the patient weans on that device.
Cfold_ETT   = cpETT.L_ETT   / nf.l_high;
Cfold_TRACH = cpTrach.L_TRACH / nf.l_high;

plot(cpTrach.cap.t, cpTrach.cap.C, '-', 'Color', c.TRACH, 'LineWidth',1.8, ...
    'DisplayName','C_{max}(t), tracheostomy (spontaneous)');
plot(cpETT.cap.t, cpETT.cap.C, '-', 'Color', c.ETT, 'LineWidth',1.8, ...
    'DisplayName','C_{max}(t), ETT (sedated)');
% HandleVisibility off: these threshold lines are annotated by their own
% Label and must not leak into the legend as 'data1'/'data2'.
yline(Cfold_TRACH, ':', 'Color', c.TRACH, 'Label','weanable on trach above here', ...
    'FontSize',7, 'HandleVisibility','off');
yline(Cfold_ETT,   ':', 'Color', c.ETT,   'Label','weanable on ETT above here', ...
    'FontSize',7, 'HandleVisibility','off');
xlabel('Days'); ylabel('C_{max}(t) (cmH_2O)');
title('Capacity trajectory vs the weaning threshold');
legend('Location','east','Box','off','FontSize',7);
ylim([0 cpTrach.p.C_max0*1.05]); grid on; set(gca,'GridAlpha',0.08);

% ---- Panel B: is the operating point weanable over time? ----
nexttile; hold on;
% For each arm, l(t) = L / C_max(t); weanable when l(t) < l_high.
lTrach = cpTrach.L_TRACH ./ cpTrach.cap.C;
lETT   = cpETT.L_ETT   ./ cpETT.cap.C;
plot(cpTrach.cap.t, lTrach, '-', 'Color', c.TRACH, 'DisplayName','tracheostomy');
plot(cpETT.cap.t,   lETT,   '-', 'Color', c.ETT,   'DisplayName','ETT (sedated)');
yline(nf.l_high, '-', 'Color', c.fold, 'LineWidth',1.4, 'Label','fold (l_{high})', 'FontSize',8, 'HandleVisibility','off');
xlabel('Days'); ylabel('operating point l(t) = L / C_{max}(t)');
title('Below the fold line = weanable');
subtitle('both axes moving at once', 'FontSize',8, 'Color', c.axisGrey);
legend('Location','northeast','Box','off','FontSize',7);
grid on; set(gca,'GridAlpha',0.08);

% ---- Panel C: window-open fraction, and the separation check ----
nexttile; hold on;
openTrach = mean(cpTrach.traj.windowOpen);
openETT   = mean(cpETT.traj.windowOpen);
openVals  = [openETT, openTrach] * 100;
b = bar(openVals, 0.6, 'FaceColor','flat', 'EdgeColor','none');
b.CData = [c.ETT; c.TRACH];
xticks([1 2]); xticklabels({'ETT (sedated)','trach (spont.)'});
ylabel('% of days weanable'); ylim([0 105]);
for i = 1:2
    text(i, openVals(i), sprintf('%.0f%%', openVals(i)), ...
        'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',9);
end
title('Time spent weanable');
subtitle(sprintf('timescale separation %.0f-%.0fx (quasi-static OK)', ...
    min(cpTrach.separation,cpETT.separation), max(cpTrach.separation,cpETT.separation)), ...
    'FontSize',8, 'Color', c.axisGrey);
grid on; set(gca,'GridAlpha',0.08);

title(tl, 'H4: the device acts on BOTH axes over time -- capacity (2b) and load (1) -- and the rescue window is their product', ...
    'FontWeight','bold','FontSize',11.5);

viz.save(fig, 'fig_vidd_dynamic_rescue');

% ---- Table & numbers ----
T = table(cpTrach.cap.t', cpTrach.cap.C', lTrach', cpETT.cap.C', lETT', ...
    'VariableNames', {'t_days','Cmax_trach','l_trach','Cmax_ETT','l_ETT'});
writetable(T, fullfile(wob.projectRoot,'results','tables','vidd_exp4_dynamic_rescue.csv'));

fprintf('\n  Timescale separation: trach %.0fx, ETT %.0fx (quasi-static coupling valid above ~5x)\n', ...
    cpTrach.separation, cpETT.separation);
fprintf('  dL_device (load axis)  = %.2f cmH2O -> rescue window %.2f cmH2O of C_max\n', ...
    cpTrach.dL_device, cpTrach.dL_device/nf.l_high);
fprintf('  capacity gap at steady state = %.1f cmH2O of C_max (trach A=%.2f vs ETT A=%.2f)\n', ...
    cpTrach.cap.Cstar_inst(end) - cpETT.cap.Cstar_inst(end), cpTrach.A, cpETT.A);
fprintf('  %% of days weanable: ETT(sedated) %.0f%%, trach(spontaneous) %.0f%%\n', 100*openETT, 100*openTrach);
fprintf('\n  CAVEAT: the capacity axis is driven by the support->A mapping, which is NOT\n');
fprintf('  calibrated. The relative size of the two axes is a consequence of that mapping.\n');
end
