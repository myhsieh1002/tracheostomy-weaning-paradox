function T = co2_exp2_requiredVE(cfgCO2, cfgM1, cfgM2)
%CO2_EXP2_REQUIREDVE H2 -- ventilation required to hold PaCO2, and its cost.
%
%   T = co2_exp2_requiredVE() holds PaCO2 at target and reports the
%   ventilation each arm demands -- and, through Model 1, what delivering it
%   costs in effort.
%
%   THIS IS THE PHYSIOLOGICALLY REAL QUESTION
%   -----------------------------------------
%   PaCO2 is regulated, so ventilation moves. Chadda 2002 measured exactly
%   this: PaCO2 and RR unchanged, V_T 330 -> 400 mL on mouth breathing. The
%   patient pays for dead space in ventilation, and ventilation is paid for
%   in work -- which is where this chains into Model 1.
%
%   The end-to-end chain is the point: CO2 kinetics set V_E, Model 1 turns
%   V_E into effort, Model 2 takes the effort as load. The right-hand panel
%   shows the whole chain in one axis.

arguments
    cfgCO2 (1,1) struct = co2.loadConfig()
    cfgM1  (1,1) struct = wob.loadConfig()
    cfgM2  (1,1) struct = lc.loadConfig()
end

c = viz.style();
arms  = {cfgCO2.arms.native, cfgCO2.arms.ett, cfgCO2.arms.trach};
label = {'not intubated', 'ETT 7.5', 'trach 8.0'};
cols  = [c.grey; c.ETT; c.TRACH];

VCO2s = linspace(150, 450, 60);

VE = nan(numel(arms), numel(VCO2s));
L  = nan(numel(arms), numel(VCO2s));
for a = 1:numel(arms)
    for k = 1:numel(VCO2s)
        cp = co2.coupling(cfgCO2, cfgM1, cfgM2, arms{a}, VCO2=VCO2s(k));
        VE(a,k) = cp.vent.V_E;
        L(a,k)  = cp.L_total;
    end
end

fig = figure('Position',[70 70 1180 430]);
tl = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

% ---- required V_E ----
nexttile; hold on;
for a = 1:numel(arms)
    plot(VCO2s, VE(a,:), '-', 'Color', cols(a,:), 'DisplayName', label{a});
end
xlabel('V̇CO_2 (mL/min)'); ylabel('Required V_E (L/min)');
title(sprintf('Ventilation to hold PaCO_2 = %g mmHg', cfgCO2.targets.PaCO2_target_mmHg));
legend('Location','northwest','Box','off');
grid on; set(gca,'GridAlpha',0.08);

% ---- the saving, absolute and relative ----
nexttile; hold on;
yyaxis left;
plot(VCO2s, VE(1,:)-VE(3,:), '-', 'Color', c.TRACH, 'DisplayName','not intubated \rightarrow trach');
plot(VCO2s, VE(2,:)-VE(3,:), '-', 'Color', c.ETT,   'DisplayName','ETT \rightarrow trach');
ylabel('\DeltaV_E saved (L/min)'); set(gca,'YColor',c.stable);
yyaxis right;
plot(VCO2s, 100*(VE(2,:)-VE(3,:))./VE(2,:), '--', 'Color', c.ETT, 'HandleVisibility','off');
ylabel('device step, % of demand'); set(gca,'YColor',c.ETT);
xlabel('V̇CO_2 (mL/min)');
title('The saving is constant; the demand is not');
subtitle('so its share falls -- dashed, right axis', 'FontSize',8, 'Color', c.axisGrey);
yyaxis left; legend('Location','northwest','Box','off');
grid on; set(gca,'GridAlpha',0.08);

% ---- the chain: CO2 -> ventilation -> load ----
nexttile; hold on;
for a = 1:numel(arms)
    plot(VCO2s, L(a,:), '-', 'Color', cols(a,:), 'DisplayName', label{a});
end
v = cfgM1.load_scale_validation.vassilakopoulos1998;
yline(v.Pi_failure_cmH2O, ':', 'Color', c.axisGrey, 'Label','P_i, weaning failure', ...
    'FontSize',8, 'HandleVisibility','off');
yline(v.Pi_success_cmH2O, ':', 'Color', c.axisGrey, 'Label','P_i, weaning success', ...
    'FontSize',8, 'HandleVisibility','off');
xlabel('V̇CO_2 (mL/min)'); ylabel('Load L = mean inspiratory P_{mus} (cmH_2O)');
title('The chain: CO_2 \rightarrow ventilation \rightarrow load');
subtitle('via Model 1; observed weaning band marked', 'FontSize',8, 'Color', c.axisGrey);
legend('Location','northwest','Box','off');
grid on; set(gca,'GridAlpha',0.08);

title(tl, 'H2: the patient pays for dead space in ventilation, and for ventilation in work', ...
    'FontWeight','bold','FontSize',12);

viz.save(fig, 'fig_co2_required_ve');

% ---- Table ----
rows = {};
for a = 1:numel(arms)
    for vco2 = cfgCO2.metabolism.VCO2_grid
        cp = co2.coupling(cfgCO2, cfgM1, cfgM2, arms{a}, VCO2=vco2);
        rows(end+1,:) = {label{a}, vco2, cp.vent.V_D*1e3, cp.vent.V_T, cp.vent.V_E, ...
            cp.vent.V_A, cp.L_total, cp.f_device}; %#ok<AGROW>
    end
end
T = cell2table(rows, 'VariableNames', ...
    {'arm','VCO2','Vd_mL','V_T_L','V_E','V_A','L_total','f_device'});
writetable(T, fullfile(wob.projectRoot,'results','tables','co2_exp2_required_ve.csv'));

fprintf('\n  Required ventilation to hold PaCO2 = %g mmHg:\n', cfgCO2.targets.PaCO2_target_mmHg);
disp(T);

% Does raising metabolism ALONE reach the observed weaning-failure load?
maxL = max(L(2,:));
if maxL >= v.Pi_failure_cmH2O
    vco2Needed = interp1(L(2,:), VCO2s, v.Pi_failure_cmH2O, 'linear');
    fprintf('  ETT load reaches the observed weaning-failure value (%.1f cmH2O) at VCO2 = %.0f mL/min\n', ...
        v.Pi_failure_cmH2O, vco2Needed);
else
    fprintf(['  ETT load peaks at %.1f cmH2O even at VCO2 = %.0f mL/min -- BELOW the observed\n' ...
             '  weaning-failure value of %.1f. Metabolic rate ALONE does not create a failing\n' ...
             '  patient in this lung (C_rs = %.2f, R_aw = %g, i.e. near-normal mechanics).\n' ...
             '  It takes the disease axes too, which is Model 1''s point: the load is a product\n' ...
             '  of mechanics AND metabolism, and neither alone gets you there.\n'], ...
        maxL, VCO2s(end), v.Pi_failure_cmH2O, cfgM1.patient.C_rs, cfgM1.patient.R_aw_native);

    % What does it take? Add the mechanical severity axis back in.
    sev = struct('C_rs', 0.02, 'R_aw_native', 15);
    Lsev = arrayfun(@(vc) co2.coupling(cfgCO2, cfgM1, cfgM2, cfgCO2.arms.ett, ...
                    VCO2=vc, PatientOverrides=sev).L_total, VCO2s);
    if max(Lsev) >= v.Pi_failure_cmH2O
        vco2Sev = interp1(Lsev, VCO2s, v.Pi_failure_cmH2O, 'linear');
        fprintf('  With severe mechanics (C_rs = %.2f, R_aw = %g) it reaches %.1f cmH2O at VCO2 = %.0f.\n', ...
            sev.C_rs, sev.R_aw_native, v.Pi_failure_cmH2O, vco2Sev);
    else
        fprintf('  Even with severe mechanics (C_rs = %.2f, R_aw = %g) it peaks at %.1f cmH2O.\n', ...
            sev.C_rs, sev.R_aw_native, max(Lsev));
    end
end
end
