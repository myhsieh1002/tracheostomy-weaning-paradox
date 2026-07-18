function T = co2_exp4_transient(cfgCO2, cfgM1)
%CO2_EXP4_TRANSIENT PaCO2 kinetics after switching ETT -> tracheostomy.
%
%   T = co2_exp4_transient() integrates the two-compartment model through the
%   moment of tracheostomy, under two behavioural assumptions.
%
%   THE TWO ASSUMPTIONS ARE THE POINT
%   ---------------------------------
%   What the transient looks like depends entirely on what the patient does
%   with the ventilation that dead space just freed up:
%
%     * FIXED VENTILATION -- V_E held where it was. PaCO2 falls, settles at a
%       new lower steady state. This is the spec's H1 read dynamically.
%     * REGULATED PaCO2 -- the patient lets V_E fall to hold PaCO2 constant.
%       Nothing happens to the gas; the benefit is taken entirely as reduced
%       work. This is what Chadda 2002 observed: PaCO2 and RR unchanged,
%       V_T adjusted.
%
%   Reality sits between, but the measured case sits at the second. Which
%   matters for how the result should be read: the dead-space benefit of a
%   tracheostomy shows up in EFFORT, not in blood gas -- so looking for it in
%   PaCO2 (as the negative studies effectively did) finds nothing, and that
%   absence is not evidence the benefit is absent.
%
%   The two compartments earn their keep here: the lung equilibrates in
%   ~35 s but the body CO2 stores take ~35 min, so the approach to the new
%   steady state is biphasic and a single-compartment model would call it
%   done long before the patient is.

arguments
    cfgCO2 (1,1) struct = co2.loadConfig()
    cfgM1  (1,1) struct = wob.loadConfig()
end

c = viz.style();
phi = co2.alveolarDeadSpaceFraction(cfgCO2, cfgM1);

tSwitch = 30;                       % min
tEnd    = 240;
tSpan   = linspace(0, tEnd, 4000);

vE = co2.requiredVentilation(cfgCO2, cfgM1, cfgCO2.arms.ett,   Phi=phi);
vT = co2.requiredVentilation(cfgCO2, cfgM1, cfgCO2.arms.trach, Phi=phi);

% Arm A -- ventilation held at the ETT value; the tube changes, V_A rises
% because the same V_E now clears more dead space.
V_A_fixedVE = @(t) ternary(t < tSwitch, vE.V_A, (vE.V_T - vT.ds.total) * vE.RR);

% Arm B -- the patient regulates PaCO2: V_A stays exactly where it was, and
% the saving is taken as less ventilation (hence less work).
V_A_regulated = @(t) vE.V_A;

outA = co2.twoCompartment(cfgCO2, tSpan, V_A_fixedVE);
outB = co2.twoCompartment(cfgCO2, tSpan, V_A_regulated);

fig = figure('Position',[70 70 1180 440]);
tl = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

% ---- PaCO2 traces ----
nexttile; hold on;
plot(outA.t, outA.PA, '-', 'Color', c.TRACH, 'DisplayName','V_E held fixed');
plot(outB.t, outB.PA, '-', 'Color', c.ETT,   'DisplayName','PaCO_2 regulated (Chadda)');
xline(tSwitch, ':', 'Color', c.fold, 'Label','tracheostomy', 'FontSize',8, 'HandleVisibility','off');
xlabel('Time (min)'); ylabel('PaCO_2 (mmHg)');
title('What happens to the blood gas');
legend('Location','northeast','Box','off');
grid on; set(gca,'GridAlpha',0.08);

% ---- arterial vs venous, showing the two timescales ----
nexttile; hold on;
plot(outA.t, outA.PA, '-',  'Color', c.TRACH, 'DisplayName','P_ACO_2 (lung, fast)');
plot(outA.t, outA.Pv, '--', 'Color', c.stable,'DisplayName','P_vCO_2 (body stores, slow)');
xline(tSwitch, ':', 'Color', c.fold, 'HandleVisibility','off');
xlabel('Time (min)'); ylabel('PCO_2 (mmHg)');
title('Why two compartments');
subtitle(sprintf('\\tau_{lung} = %.2f min, \\tau_{body} = %.0f min (ratio %.0f)', ...
    outA.tau_lung, outA.tau_body, outA.tau_body/outA.tau_lung), ...
    'FontSize',8, 'Color', c.axisGrey);
legend('Location','northeast','Box','off');
grid on; set(gca,'GridAlpha',0.08);

% ---- where the benefit actually lands ----
nexttile; hold on;
cats = categorical({'\DeltaPaCO_2 (mmHg)','\DeltaV_E (L/min)'});
cats = reordercats(cats, {'\DeltaPaCO_2 (mmHg)','\DeltaV_E (L/min)'});
dA = [outA.PA(end) - outA.PA(1), 0];
dB = [outB.PA(end) - outB.PA(1), vT.V_E - vE.V_E];
b = bar(cats, [dA; dB]', 'EdgeColor','none');
b(1).FaceColor = c.TRACH; b(1).DisplayName = 'V_E held fixed';
b(2).FaceColor = c.ETT;   b(2).DisplayName = 'PaCO_2 regulated';
yline(0, '-', 'Color', c.stable, 'HandleVisibility','off');
ylabel('change after tracheostomy');
title('Where the dead-space benefit lands');
subtitle('in gas, or in ventilation — never in both', 'FontSize',8, 'Color', c.axisGrey);
legend('Location','southwest','Box','off');
grid on; set(gca,'GridAlpha',0.08);

title(tl, 'The dead-space benefit of a tracheostomy shows up in effort, not in blood gas — which is why measuring PaCO_2 finds nothing', ...
    'FontWeight','bold','FontSize',11.5);

viz.save(fig, 'fig_co2_transient');

T = table(outA.t', outA.PA', outA.Pv', outB.PA', outB.Pv', ...
    'VariableNames', {'t_min','PaCO2_fixedVE','PvCO2_fixedVE','PaCO2_regulated','PvCO2_regulated'});
writetable(T, fullfile(wob.projectRoot,'results','tables','co2_exp4_transient.csv'));

fprintf('\n  Timescales: tau_lung = %.2f min, tau_body = %.1f min (ratio %.0f)\n', ...
    outA.tau_lung, outA.tau_body, outA.tau_body/outA.tau_lung);
fprintf('  Fixed-V_E arm:   PaCO2 %.1f -> %.1f mmHg (%+.2f)\n', ...
    outA.PA(1), outA.PA(end), outA.PA(end)-outA.PA(1));
fprintf('  Regulated arm:   PaCO2 %.1f -> %.1f mmHg (%+.2f), V_E %.2f -> %.2f L/min (%+.2f)\n', ...
    outB.PA(1), outB.PA(end), outB.PA(end)-outB.PA(1), vE.V_E, vT.V_E, vT.V_E-vE.V_E);
fprintf('\n  Chadda 2002 observed the regulated case: PaCO2 and RR unchanged, V_T adjusted.\n');
fprintf('  So a study looking for the tracheostomy''s dead-space benefit in PaCO2 will find\n');
fprintf('  nothing -- not because the benefit is absent, but because it is taken as work.\n');

% Half-time to the new steady state, to show the body store dominates.
if abs(outA.PA(end) - outA.PA(1)) > 1e-3
    target = outA.PA(1) + 0.5*(outA.PA(end) - outA.PA(1));
    post = outA.t > tSwitch;
    tHalf = interp1(outA.PA(post), outA.t(post), target, 'linear', NaN) - tSwitch;
    fprintf('\n  Half-time to the new PaCO2 steady state: %.1f min -- set by the body stores,\n', tHalf);
    fprintf('  not the lung. A single-compartment model would report ~%.1f min and be wrong.\n', outA.tau_lung*0.69);
end
end

function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end
