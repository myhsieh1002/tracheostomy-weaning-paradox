function T = vidd_exp3_catabolic(cfgV, cfgM1)
%VIDD_EXP3_CATABOLIC H3, reframed by the data: offset vs drain.
%
%   T = vidd_exp3_catabolic() contrasts the two ways sepsis can enter the
%   model, because the literature and the spec disagree about which is real.
%
%   THE SPEC'S H3 IS CONTRADICTED, AND THE FIGURE SHOWS WHY
%   ------------------------------------------------------
%   The spec's H3: disease catabolism is a persistent, activity-independent
%   drain (d_disease*C) that no ventilation strategy can rescue, so capacity
%   collapses regardless of how well the patient is ventilated.
%
%   The best human data say the opposite:
%     * Demoule 2013 -- sepsis is an OFFSET present at admission (-3.74
%       cmH2O), with no measurable Day1->Day3 decline. An initial condition,
%       not a slope.
%     * Lecronier 2022 -- septic patients' diaphragm force ROSE 19% over 4
%       days while ventilated, vs -7% in non-septic patients. The title is
%       "Severe but REVERSIBLE".
%
%   So the panel runs BOTH:
%     d_mode = 'offset' (data)  -- sepsis lowers C(0); the patient then
%       RECOVERS along the same trajectory a non-septic patient follows,
%       just from a lower start. Strategy still matters, and sepsis does not
%       abolish rescue.
%     d_mode = 'rate' (spec)    -- a true drain; capacity collapses and,
%       above a threshold d_disease, no activity can hold it.
%
%   The contrast IS the result: which mechanism the model uses changes the
%   clinical conclusion, and only one of them is supported.

arguments
    cfgV  (1,1) struct = vidd.loadConfig()
    cfgM1 (1,1) struct = wob.loadConfig()
end

c = viz.style();
p = vidd.calibrate(cfgV);
days = cfgV.sim.days;
A_trach = vidd.supportToActivity(cfgV.strategy.scenarios.trach_spontaneous.support_level, cfgV);
A_sed   = vidd.supportToActivity(cfgV.strategy.scenarios.deep_sedation_controlled.support_level, cfgV);

fig = figure('Position',[60 60 1240 440]);
tl = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

% ---- Panel A: offset (data) -- septic starts low, recovers ----
nexttile; hold on;
for sept = [false true]
    rT = vidd.simulateCapacity(cfgV, p, @(t) A_trach, IsSeptic=sept, Days=days);
    rS = vidd.simulateCapacity(cfgV, p, @(t) A_sed,   IsSeptic=sept, Days=days);
    style = '-'; if sept, style = '--'; end
    plot(rT.t, rT.C, style, 'Color', c.TRACH, 'DisplayName', sprintf('spontaneous, septic=%d', sept));
    plot(rS.t, rS.C, style, 'Color', c.ETT,   'DisplayName', sprintf('sedated, septic=%d', sept));
end
xlabel('Days'); ylabel('C_{max}(t) (cmH_2O)');
title('d\_mode = offset (data)');
subtitle('sepsis = lower start, then recovery', 'FontSize',8, 'Color', c.rescue);
legend('Location','east','Box','off','FontSize',7);
ylim([0 p.C_max0*1.05]); grid on; set(gca,'GridAlpha',0.08);

% ---- Panel B: rate (spec) -- collapse regardless of strategy ----
nexttile; hold on;
dgrid = cfgV.disease.d_disease_grid;
for k = 1:numel(dgrid)
    q = p; q.d_mode = "rate"; q.d_disease = dgrid(k);
    rT = vidd.simulateCapacity(cfgV, q, @(t) A_trach, Days=days);
    shade = 0.2 + 0.8*(k-1)/max(numel(dgrid)-1,1);
    plot(rT.t, rT.C, '-', 'Color', c.ETT*shade + [1 1 1]*(1-shade), ...
        'DisplayName', sprintf('d = %.2f', dgrid(k)));
end
xlabel('Days'); ylabel('C_{max}(t) (cmH_2O)');
title('d\_mode = rate (spec)');
subtitle('a true drain: strategy cannot hold it', 'FontSize',8, 'Color', c.axisGrey);
legend('Location','northeast','Box','off','FontSize',7);
ylim([0 p.C_max0*1.05]); grid on; set(gca,'GridAlpha',0.08);

% ---- Panel C: strategy leverage vs catabolism (H3 quantified) ----
nexttile; hold on;
dsweep = linspace(0, 0.3, 40);
absLev = zeros(size(dsweep)); relLev = zeros(size(dsweep));
for k = 1:numel(dsweep)
    q = p; q.d_mode = "rate"; q.d_disease = dsweep(k);
    l = vidd.strategyLeverage(q);
    absLev(k) = l.absolute; relLev(k) = l.relative;
end
yyaxis left;
plot(dsweep, absLev, '-', 'Color', c.stable, 'LineWidth',1.8, 'DisplayName','absolute (cmH_2O)');
ylabel('strategy leverage, absolute (cmH_2O)'); set(gca,'YColor',c.stable);
yyaxis right;
plot(dsweep, relLev, '--', 'Color', c.grey, 'DisplayName','relative (ratio)');
ylabel('relative (C_{best}/C_{worst})'); set(gca,'YColor',c.grey);
xlabel('d\_disease (1/day)');
title('How much can strategy still buy?');
subtitle('absolute leverage collapses; relative does not', 'FontSize',8, 'Color', c.axisGrey);
yyaxis left; legend('Location','east','Box','off','FontSize',7);
grid on; set(gca,'GridAlpha',0.08);

title(tl, 'H3 reframed: the data make sepsis a reversible starting offset, not the drain the spec assumed', ...
    'FontWeight','bold','FontSize',11.5);

viz.save(fig, 'fig_vidd_catabolic');

% ---- Table ----
rows = {};
for sept = [false true]
    for str = ["trach_spontaneous","deep_sedation_controlled"]
        A = vidd.supportToActivity(cfgV.strategy.scenarios.(str).support_level, cfgV);
        r = vidd.simulateCapacity(cfgV, p, @(t) A, IsSeptic=sept, Days=days);
        rows(end+1,:) = {'offset', char(str), sept, 0, r.C0, r.C(end)}; %#ok<AGROW>
    end
end
for d = cfgV.disease.d_disease_grid
    q = p; q.d_mode = "rate"; q.d_disease = d;
    r = vidd.simulateCapacity(cfgV, q, @(t) A_trach, Days=days);
    rows(end+1,:) = {'rate', 'trach_spontaneous', false, d, r.C0, r.C(end)}; %#ok<AGROW>
end
T = cell2table(rows, 'VariableNames', {'d_mode','scenario','septic','d_disease','C0','C_end'});
writetable(T, fullfile(wob.projectRoot,'results','tables','vidd_exp3_catabolic.csv'));

fprintf('\n  H3 VERDICT: the spec''s persistent-drain mechanism (d_mode=rate) is what makes\n');
fprintf('  capacity uncatchable. The data (d_mode=offset) make sepsis a reversible offset:\n');
rOff = vidd.simulateCapacity(cfgV, p, @(t) A_trach, IsSeptic=true, Days=days);
fprintf('    septic on spontaneous breathing: C0=%.1f -> C(%dd)=%.1f (RECOVERS)\n', rOff.C0, days, rOff.C(end));
fprintf('  Lecronier 2022 observed exactly this: septic diaphragm force +19%% over 4 d.\n');
end
