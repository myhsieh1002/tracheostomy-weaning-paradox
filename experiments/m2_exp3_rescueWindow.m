function T = m2_exp3_rescueWindow(cfg1, cfg2)
%M2_EXP3_RESCUEWINDOW H3 -- where a tracheostomy can flip the outcome. THESIS FIGURE.
%
%   T = m2_exp3_rescueWindow() maps the region in which switching from an
%   ETT to a tracheostomy moves the operating point across the fold, with
%   loads taken end-to-end from Model 1 (never hard-coded).
%
%   THREE FINDINGS THAT DEPART FROM THE BUILD SPEC
%   ----------------------------------------------
%   1. The window does not SHRINK with capacity -- it MOVES. Because
%
%          window width in C_max = dL_device / l_high
%
%      and dL_device is near-constant across severity, the band keeps its
%      width and slides upward as disease raises the load. The spec
%      predicted a window collapsing to nothing at low C_max; that is not
%      what the model does.
%
%   2. Metabolic rate, which the spec did not treat as a severity axis at
%      all, dominates where the window sits. At a normal V_A of 4.2 L/min
%      the model cannot even reach the loads measured in failing weaning
%      patients (max 13.6 vs Vassilakopoulos's 19.5 cmH2O), and the window
%      falls below the observed MIP range entirely -- a spurious "the device
%      never matters" result that is an artefact of assuming normal
%      metabolism. Only at V_A ~ 8.4 (VCO2 ~ 400 mL/min, severe catabolism)
%      does the window overlap the population the trials actually enrol.
%
%   3. The window is narrow wherever it sits: ~2-5 cmH2O of C_max against a
%      clinical range of 20-70, i.e. under 10% of the capacity axis.
%
%   THE PARADOX, QUANTIFIED
%   -----------------------
%   For a tracheostomy to decide a patient's outcome, that patient's
%   capacity must land inside a band a few cmH2O wide whose position depends
%   on their compliance, airway resistance AND metabolic rate -- none of
%   which is known precisely at the bedside. An unselected trial puts few
%   patients in the window and averages to approximately nothing. That is
%   the weaning-survival paradox, in numbers.

arguments
    cfg1 (1,1) struct = wob.loadConfig()
    cfg2 (1,1) struct = lc.loadConfig()
end

c = viz.style();
p = lc.calibrate(cfg2);
nf = lc.normalizedFolds(p);

devETT   = cfg2.load_from_model1.devices{1};
devTRACH = cfg2.load_from_model1.devices{2};

CmaxAxis = linspace(10, 80, 240);
Cvals    = linspace(0.02, 0.06, 110);
R_fixed  = 10;
VAlevels = cfg1.disease_grid.target_VA_L_min(:)';
VCO2     = round(VAlevels / 4.2 * 200);

vald  = cfg1.load_scale_validation.vassilakopoulos1998;
MIP   = [42.3, 53.8];   % observed weaning-failure / success capacity

cmap = [c.TRACH_light; c.rescue; c.ETT_light];   % both_wean | rescued | both_fail

fig = figure('Position',[50 50 1320 460]);
tl = tiledlayout(1, numel(VAlevels)+1, 'TileSpacing','compact','Padding','compact');

rows = {};
for v = 1:numel(VAlevels)
    cfgV = cfg1;
    cfgV.patient.target_VA_L_min = VAlevels(v);

    L_E = zeros(size(Cvals)); L_T = zeros(size(Cvals));
    for j = 1:numel(Cvals)
        ov = struct('C_rs', Cvals(j), 'R_aw_native', R_fixed);
        L_E(j) = lc.coupling(cfgV, cfg2, devETT,   ov).L_total;
        L_T(j) = lc.coupling(cfgV, cfg2, devTRACH, ov).L_total;
    end

    M = classifyGrid(L_E, L_T, CmaxAxis, nf.l_high);

    nexttile; hold on;
    imagesc(Cvals, CmaxAxis, M); set(gca,'YDir','normal');
    colormap(gca, cmap); clim([1 3]); axis tight;

    % Analytic window boundaries
    plot(Cvals, L_E ./ nf.l_high, '-', 'Color','w', 'LineWidth',1.4);
    plot(Cvals, L_T ./ nf.l_high, '-', 'Color','w', 'LineWidth',1.4);

    % The capacity band the weaning literature actually reports
    patch([Cvals(1) Cvals(end) Cvals(end) Cvals(1)], [MIP(1) MIP(1) MIP(2) MIP(2)], ...
        'k', 'FaceAlpha',0.10, 'EdgeColor','k', 'LineStyle',':', 'LineWidth',1);
    text(Cvals(2), mean(MIP), ' observed MIP range', 'FontSize',8, 'Color','k');

    xlabel('C_{rs} (L/cmH_2O)  \rightarrow  healthier');
    if v == 1, ylabel('C_{max} (cmH_2O)'); end
    title(sprintf('V_A = %.1f L/min  (VCO_2 \\approx %d)', VAlevels(v), VCO2(v)));

    overlaps = any(L_E./nf.l_high > MIP(1) & L_T./nf.l_high < MIP(2));
    subtitle(ternary(overlaps, 'window reaches the observed population', ...
                               'window sits BELOW the observed population'), ...
             'FontSize',9, 'Color', ternary(overlaps, c.rescue, c.grey));

    for j = 1:numel(Cvals)
        rows(end+1,:) = {VAlevels(v), VCO2(v), Cvals(j), R_fixed, L_E(j), L_T(j), ...
            L_E(j)-L_T(j), L_T(j)/nf.l_high, L_E(j)/nf.l_high, (L_E(j)-L_T(j))/nf.l_high, overlaps}; %#ok<AGROW>
    end
end

% ---- Final panel: window width vs dL_device, the analytic relation ----
nexttile; hold on;
dLAxis = linspace(0, 4, 100);
plot(dLAxis, dLAxis / nf.l_high, '-', 'Color', c.stable, 'LineWidth', 2, ...
    'DisplayName', sprintf('width = \\DeltaL / %.3f', nf.l_high));

T_tmp = cell2table(rows, 'VariableNames', ...
    {'target_VA','VCO2','C_rs','R_aw','L_ETT','L_TRACH','dL_device', ...
     'window_Cmax_low','window_Cmax_high','window_width','overlaps_observed_MIP'});

for v = 1:numel(VAlevels)
    sel = T_tmp.target_VA == VAlevels(v);
    scatter(T_tmp.dL_device(sel), T_tmp.window_width(sel), 12, ...
        'MarkerFaceColor', c.TRACH * (0.4 + 0.3*v), 'MarkerEdgeColor','none', ...
        'DisplayName', sprintf('V_A = %.1f', VAlevels(v)));
end
yline(diff([20 70]), '--', 'Color', c.grey, ...
    'Label','full clinical C_{max} range', 'FontSize',8, 'HandleVisibility','off');

xlabel('\DeltaL_{device} = L_{ETT} - L_{TRACH}  (cmH_2O)');
ylabel('Rescue window width in C_{max} (cmH_2O)');
title('The window is only as wide as the load removed');
legend('Location','northwest','Box','off');
grid on; set(gca,'GridAlpha',0.08);

h(1) = patch(NaN,NaN,cmap(1,:),'DisplayName','both wean (device not limiting)');
h(2) = patch(NaN,NaN,cmap(2,:),'DisplayName','RESCUED by tracheostomy');
h(3) = patch(NaN,NaN,cmap(3,:),'DisplayName','both fail (capacity too low)');
lg = legend(h, 'Orientation','horizontal','Box','off');
lg.Layout.Tile = 'south';

title(tl, 'H3: the tracheostomy decides the outcome only inside a narrow, disease-dependent band of capacity', ...
    'FontWeight','bold','FontSize',13);

viz.save(fig, 'fig_rescue_window');

T = T_tmp;
writetable(T, fullfile(wob.projectRoot,'results','tables','m2_exp3_rescue.csv'));

% ================= Reported numbers =================
fprintf('\n  l_high = %.4f  (normalised fold; window width = dL_device / l_high)\n\n', nf.l_high);
fprintf('  %6s %6s %10s %8s %8s %8s %18s %s\n', 'V_A','VCO2','C_rs','L_ETT','L_TRACH','dL','window C_max','reaches obs?');
for v = 1:numel(VAlevels)
    for cr = [0.02 0.04 0.06]
        [~, j] = min(abs(Cvals - cr));
        sel = find(T.target_VA == VAlevels(v) & T.C_rs == Cvals(j), 1);
        fprintf('  %6.1f %6d %10.3f %8.2f %8.2f %8.2f    (%5.1f, %5.1f)  %s\n', ...
            T.target_VA(sel), T.VCO2(sel), T.C_rs(sel), T.L_ETT(sel), T.L_TRACH(sel), ...
            T.dL_device(sel), T.window_Cmax_low(sel), T.window_Cmax_high(sel), ...
            ternary(T.overlaps_observed_MIP(sel),'YES',''));
    end
end

fprintf('\n  Load-scale check against Vassilakopoulos (Pi = %.1f fail / %.1f success):\n', ...
    vald.Pi_failure_cmH2O, vald.Pi_success_cmH2O);
fprintf('    max L_ETT reachable: %.2f cmH2O at V_A=%.1f -> %s\n', max(T.L_ETT), max(VAlevels), ...
    ternary(max(T.L_ETT) >= vald.Pi_success_cmH2O, 'reaches the observed range', 'BELOW the observed range'));
fprintf('    mean window width %.2f cmH2O = %.0f%% of a 20-70 capacity range\n', ...
    mean(T.window_width), 100*mean(T.window_width)/50);
end

% ================= helpers =================

function M = classifyGrid(L_E, L_T, CmaxAxis, l_high)
%CLASSIFYGRID Outcome code over a (C_max x severity) grid, vectorised.
%
%   Codes: 1 = both wean, 2 = rescued by tracheostomy, 3 = both fail.
%
%   A sustainable branch exists iff L < l_high*C_max. Once l_high is known
%   the classification is elementwise arithmetic, so the grid needs no loop
%   and no per-cell branch solve. l_high is computed ONCE by the caller:
%   by the scale invariance it depends only on (beta/alpha, s, u_crit), so
%   recomputing it per cell would repeat an identical ~0.08 s solve tens of
%   thousands of times.
%
%   lc.rescueOutcome applies the identical rule one patient at a time;
%   tests/tRescueApi.m checks the two agree.

L_fold = l_high * CmaxAxis(:);
ettWeans   = L_fold > L_E(:)';
trachWeans = L_fold > L_T(:)';

M = zeros(numel(CmaxAxis), numel(L_E));
M(ettWeans  &  trachWeans) = 1;
M(~ettWeans &  trachWeans) = 2;
M(~ettWeans & ~trachWeans) = 3;
M(ettWeans  & ~trachWeans) = NaN;   % impossible while L_T < L_E
end

function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end
