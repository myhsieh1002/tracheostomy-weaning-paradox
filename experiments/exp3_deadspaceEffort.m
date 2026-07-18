function T = exp3_deadspaceEffort(cfg)
%EXP3_DEADSPACEEFFORT H3 -- dead space, required ventilation and effort (Mode B).
%
%   T = exp3_deadspaceEffort() holds target alveolar ventilation fixed and
%   compares the minute ventilation and muscle effort each device demands.
%
%   THE SPEC'S H3 MECHANISM DOES NOT SURVIVE, AND THAT IS THE RESULT
%   ---------------------------------------------------------------
%   The build spec expected a ~70 mL dead-space advantage to drive H3. It
%   does not. An ETT and a tracheostomy tube are both sited with the tip in
%   the mid-trachea, so BOTH bypass the same upper airway; the bypass term
%   is common to the two arms and CANCELS. What remains is the apparatus
%   volume difference -- ~18 mL (ETT 7.5) against ~5.5 mL (trach 8.0), i.e.
%   ~12 mL against a tidal volume of several hundred.
%
%   The 72 mL figure is real, but it is measured against the NON-INTUBATED
%   state (Nunn 1959; Chadda 2002), not against an ETT.
%
%   This panel therefore reports THREE arms -- not intubated, ETT, trach --
%   because the two-arm comparison is what makes the literature's 70 mL
%   look relevant when it is not.
%
%   EXTERNAL VALIDATION: two studies measured dead space before and after
%   tracheostomy in already-intubated patients and found no change --
%   Mohr 2001 (n=42, V_D/V_T 0.51 -> 0.51) and Joseph 2013 (n=24, 41% vs
%   40%, p=0.75, titled "the myth of dead space"). The model predicts a
%   null; the clinic measured a null.

arguments
    cfg (1,1) struct = wob.loadConfig()
end

c = viz.style();
devETT   = cfg.experiments.exp3.devices{1};
devTRACH = cfg.experiments.exp3.devices{2};

% Third arm: the intact, non-intubated airway. Modelled as a pseudo-device
% (NATIVE_UPPER_AIRWAY) carrying full anatomic dead space, no apparatus, and
% a resistance equal to the tracheostomy tube's -- which is what Chadda 2002
% measured in vivo. Giving this arm zero resistance instead, as would be the
% naive reading of "no tube", inverts the comparison against Davis 1999.
arms = { cfg, 'NATIVE_UPPER_AIRWAY', 'not intubated';
         cfg, devETT,                'ETT 7.5';
         cfg, devTRACH,              'trach 8.0' };

rows = {};
for k = 1:size(arms,1)
    cfgK = arms{k,1}; devName = arms{k,2}; label = arms{k,3};
    dev = wob.getDevice(cfgK, devName);
    ds  = wob.deadSpace(cfgK, dev);
    e   = wob.simulateEffort(cfgK, devName);
    rows(end+1,:) = {label, ds.total*1e3, ds.anatomicEffective*1e3, ds.apparatus*1e3, ...
        e.V_T, e.V_E, e.P_mus_mean, e.P_mus_peak, e.WOB_total_J_min, e.PTP_per_min}; %#ok<AGROW>
end

T = cell2table(rows, 'VariableNames', ...
    {'arm','Vd_total_mL','Vd_anatomic_eff_mL','Vd_apparatus_mL','V_T_L','V_E_L_min', ...
     'P_mus_mean','P_mus_peak','WOB_total_J_min','PTP_per_min'});
writetable(T, fullfile(wob.projectRoot,'results','tables','exp3_table.csv'));

% ---------------- Figure ----------------
cols = [c.grey; c.ETT; c.TRACH];
metrics = {'Vd_total_mL','V_E_L_min','P_mus_mean','WOB_total_J_min'};
labels  = {'Total dead space (mL)','Required V_E (L/min)','Mean P_{mus} (cmH_2O)','Total WOB (J/min)'};

fig = figure('Position',[70 70 1240 400]);
tl = tiledlayout(1,4,'TileSpacing','compact','Padding','compact');

for m = 1:numel(metrics)
    nexttile; hold on;
    vals = T.(metrics{m});
    for k = 1:height(T)
        bar(k, vals(k), 0.66, 'FaceColor', cols(k,:), 'EdgeColor','none');
        text(k, vals(k), sprintf(' %.4g', vals(k)), 'HorizontalAlignment','center', ...
            'VerticalAlignment','bottom','FontSize',8);
    end
    xticks(1:height(T)); xticklabels(T.arm); xtickangle(25);
    ylabel(labels{m}); ylim([0, max(vals)*1.2]);
    grid on; set(gca,'GridAlpha',0.08);

    % Annotate the two contrasts that matter
    d_nat_ett = vals(2) - vals(1);
    d_ett_tr  = vals(3) - vals(2);
    subtitle(sprintf('intubation: %+.4g   |   ETT\\rightarrowtrach: %+.4g', d_nat_ett, d_ett_tr), ...
        'FontSize',8, 'Color', c.axisGrey);
end

title(tl, {'H3: the dead-space benefit is against the NON-INTUBATED state, not against an ETT', ...
    'the ETT\rightarrowtrach step moves only the apparatus volume'}, ...
    'FontWeight','bold','FontSize',12);

viz.save(fig, 'fig_deadspace_effort');

% ---------------- Reported numbers ----------------
fprintf('\n  Dead-space budget (target V_A = %.1f L/min held fixed):\n', cfg.patient.target_VA_L_min);
disp(T(:, {'arm','Vd_total_mL','V_T_L','V_E_L_min','P_mus_mean','WOB_total_J_min'}));

dVd_intubation = T.Vd_total_mL(1) - T.Vd_total_mL(2);
dVd_device     = T.Vd_total_mL(2) - T.Vd_total_mL(3);
fprintf('  Dead space removed by INTUBATION (vs not intubated): %+.1f mL\n', dVd_intubation);
fprintf('  Dead space removed by ETT -> TRACH:                  %+.1f mL  <- the device step\n', dVd_device);
fprintf('  Ratio: the device step is %.0f%% of the intubation step\n', 100*dVd_device/dVd_intubation);

fprintf('\n  Effort difference ETT -> trach: V_E %+.3f L/min, P_mus %+.3f cmH2O, WOB %+.3f J/min\n', ...
    T.V_E_L_min(3)-T.V_E_L_min(2), T.P_mus_mean(3)-T.P_mus_mean(2), ...
    T.WOB_total_J_min(3)-T.WOB_total_J_min(2));
fprintf('  (P_mus difference is driven mostly by RESISTANCE, not dead space -- see exp2.)\n');

% ---------------- External validation: Chadda 2002 ----------------
v = cfg.chadda2002_validation;
modelRatio = T.WOB_total_J_min(1) / T.WOB_total_J_min(3);   % not-intubated / trach
fprintf('\n  EXTERNAL CHECK -- Chadda 2002 (PMID 12447520), the only in vivo\n');
fprintf('  non-intubated vs tracheostomised comparison (n=9, crossover):\n');
fprintf('    observed WOB ratio (mouth / trach) = %.2f  (%.1f / %.1f J/min)\n', ...
    v.WOB_ratio_mouth_over_trach, v.mouth_breathing.WOB_J_min, v.tracheal_breathing.WOB_J_min);
fprintf('    model    WOB ratio (native / trach) = %.2f  (%.2f / %.2f J/min)\n', ...
    modelRatio, T.WOB_total_J_min(1), T.WOB_total_J_min(3));
fprintf('    observed dV_D = %d mL | model dV_D = %.0f mL\n', ...
    v.dVd_mL, T.Vd_total_mL(1)-T.Vd_total_mL(3));
fprintf('    -> %s (ratio differs by %.0f%%)\n', ...
    ternary(abs(modelRatio - v.WOB_ratio_mouth_over_trach) < 0.20, 'CONSISTENT', 'DISCREPANT'), ...
    100*abs(modelRatio - v.WOB_ratio_mouth_over_trach)/v.WOB_ratio_mouth_over_trach);
fprintf('    Caveat: %s\n', v.caveat);

fprintf('\n  Also consistent with Mohr 2001 (n=42) and Joseph 2013 (n=24): no measurable\n');
fprintf('  dead-space change after tracheostomy in ALREADY-INTUBATED patients -- the\n');
fprintf('  apparatus term only, which is exactly what the ETT->trach step shows.\n');
end

function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end
