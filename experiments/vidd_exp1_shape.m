function T = vidd_exp1_shape(cfgV, cfgM1)
%VIDD_EXP1_SHAPE H1 -- equilibrium capacity vs activity, both shapes.
%
%   T = vidd_exp1_shape() plots C*(A) under the two degradation shapes and
%   shows what each is and is not supported by.
%
%   THE HONEST VERSION OF H1
%   ------------------------
%   The spec's H1 is an inverted U in capacity vs activity, from a U-shaped
%   degradation g(A). The literature splits on this:
%
%     * MONOTONIC (default) -- calibrated to Zambon 2016, the only study
%       stratifying atrophy RATE by ventilatory support. It is monotone:
%       best at full spontaneous activity, no upturn. C*(A) rises with A.
%
%     * U-SHAPE (the spec) -- supported only in a CLINICAL OUTCOME
%       (Goligher 2018: thickening fraction 15-30% gives the shortest
%       ventilation), never in a degradation rate.
%
%   Both are drawn. The point of the panel is that the shape the model uses
%   is a CHOICE between a rate-calibrated monotone form and an
%   outcome-inspired U, and the reader should see both rather than be handed
%   one.

arguments
    cfgV  (1,1) struct = vidd.loadConfig()
    cfgM1 (1,1) struct = wob.loadConfig()
end

c = viz.style();
p = vidd.calibrate(cfgV);

A = linspace(0, 1, 400);
pMono = p; pMono.g_mode = "monotonic";
pU    = p; pU.g_mode    = "ushape";

Cmono = vidd.equilibriumCapacity(A, pMono);
Cu    = vidd.equilibriumCapacity(A, pU);
gMono = vidd.degradation(A, pMono);
gU    = vidd.degradation(A, pU);

fig = figure('Position',[70 70 1180 440]);
tl = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

% ---- degradation g(A) ----
nexttile; hold on;
plot(A, gMono, '-', 'Color', c.TRACH, 'DisplayName','monotonic (Zambon rate)');
plot(A, gU,    '-', 'Color', c.ETT,   'DisplayName','U-shape (spec)');
% Zambon's measured points (normalised to CMV)
z = cfgV.validation.zambon2016;
useNeg = z.dTdi_pct_per_day < 0;
A_z = z.assigned_A(useNeg);
g_z = z.dTdi_pct_per_day(useNeg) / z.dTdi_pct_per_day(1);
scatter(A_z, g_z, 45, c.stable, 'filled', 'MarkerEdgeColor','w', 'DisplayName','Zambon 2016 (measured)');
xline(pU.A_star, ':', 'Color', c.axisGrey, 'Label','A^* (U only)', 'FontSize',8, 'HandleVisibility','off');
xlabel('Diaphragm activity A'); ylabel('degradation g(A)');
title('Degradation rate vs activity');
legend('Location','north','Box','off');
grid on; set(gca,'GridAlpha',0.08);

% ---- equilibrium capacity C*(A) ----
nexttile; hold on;
plot(A, Cmono, '-', 'Color', c.TRACH, 'DisplayName','monotonic');
plot(A, Cu,    '-', 'Color', c.ETT,   'DisplayName','U-shape');
yline(p.C_max0, ':', 'Color', c.axisGrey, 'Label','C_{max0}', 'FontSize',8, 'HandleVisibility','off');
% Model 2's rescue band for reference
pM2 = lc.calibrate(lc.loadConfig());
% (no fold line here -- C_max is the axis, not the load; shown in exp4)
xlabel('Diaphragm activity A'); ylabel('Equilibrium capacity C^* (cmH_2O)');
title('Capacity vs activity');
subtitle('monotonic rises to the ceiling; U peaks then falls', 'FontSize',8, 'Color', c.axisGrey);
legend('Location','southeast','Box','off');
grid on; set(gca,'GridAlpha',0.08);

% ---- the h0=1 test: shape is g's, not h's ----
nexttile; hold on;
for h0 = [p.h0, 1.0]
    q = pU; q.h0 = h0;
    plot(A, vidd.equilibriumCapacity(A, q), '-', 'LineWidth', 1.4 + (h0==1)*0.4, ...
        'Color', c.ETT*(0.4 + 0.6*(h0<1)) + [1 1 1]*0.6*(h0==1)*0, ...
        'DisplayName', sprintf('U-shape, h0 = %.2f', h0));
end
xlabel('Diaphragm activity A'); ylabel('C^* (cmH_2O)');
title('The U is g''s doing, not h''s');
subtitle('the peak survives h0 = 1 (h constant)', 'FontSize',8, 'Color', c.axisGrey);
legend('Location','south','Box','off');
grid on; set(gca,'GridAlpha',0.08);

title(tl, 'H1: capacity depends on activity -- but whether it is monotone or U-shaped is a choice the data underdetermine', ...
    'FontWeight','bold','FontSize',11.5);

viz.save(fig, 'fig_vidd_shape');

% ---- Table & verdict ----
levM = vidd.strategyLeverage(pMono);
levU = vidd.strategyLeverage(pU);
T = table(A', Cmono', Cu', gMono', gU', 'VariableNames', {'A','Cstar_mono','Cstar_ushape','g_mono','g_ushape'});
writetable(T, fullfile(wob.projectRoot,'results','tables','vidd_exp1_shape.csv'));

fprintf('\n  Monotonic: peak C*=%.1f at A=%.2f, worst %.1f at A=%.2f (rises to ceiling)\n', ...
    levM.C_best, levM.A_best, levM.C_worst, levM.A_worst);
fprintf('  U-shape:   peak C*=%.1f at A=%.2f, worst %.1f at A=%.2f\n', ...
    levU.C_best, levU.A_best, levU.C_worst, levU.A_worst);
fprintf('  Under the spec''s A definition (A=1 = fully spontaneous), the DATA (Zambon)\n');
fprintf('  put the optimum at A=1, i.e. monotonic. The U''s injury arm needs A>1.\n');
end
