function T = m2_exp1_bifurcation1p(cfg2)
%M2_EXP1_BIFURCATION1P H1 -- bistability and hysteresis in load.
%
%   T = m2_exp1_bifurcation1p() draws the one-parameter bifurcation diagram
%   C* vs L, marks both saddle-node folds, and demonstrates the hysteresis
%   loop by integrating the ODE up and down through the bistable window.
%
%   The branch itself comes from the closed-form inverse
%   (lc.equilibriumBranch); the folds are its stationary points. The
%   forward root-finder (lc.fixedPoints) and a direct ode45 integration are
%   overlaid as independent checks that the closed form is right.

arguments
    cfg2 (1,1) struct = lc.loadConfig()
end

c = viz.style();
p = lc.calibrate(cfg2);
C_max = 50;

br = lc.equilibriumBranch(C_max, p);

fig = figure('Position',[80 80 1120 440]);
tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

% ---------------- Panel 1: bifurcation diagram ----------------
nexttile; hold on;

% Split the curve into stable/unstable runs so the dashes break correctly.
plotRuns(br.L, br.C, br.stable,  c.stable,   '-',  2.0, 'stable');
plotRuns(br.L, br.C, ~br.stable, c.unstable, '--', 1.4, 'unstable (separatrix)');

for k = 1:numel(br.folds)
    plot(br.folds(k).L, br.folds(k).C, 'o', 'MarkerSize',8, ...
        'MarkerFaceColor',c.fold, 'MarkerEdgeColor','w', 'LineWidth',1, ...
        'HandleVisibility', ternary(k==1,'on','off'), 'DisplayName','saddle-node fold');
end

% Bistable window shading
yl = [0, C_max*1.05];
patch([br.L_fold_low br.L_fold_high br.L_fold_high br.L_fold_low], ...
      [yl(1) yl(1) yl(2) yl(2)], c.fold, 'FaceAlpha',0.06, 'EdgeColor','none', ...
      'HandleVisibility','off');

% Independent check: forward root-finding at sampled loads
Lcheck = linspace(0.2, br.L_fold_high*1.4, 60);
for L = Lcheck
    fp = lc.fixedPoints(L, C_max, p);
    scatter(repmat(L,1,numel(fp.C)), fp.C, 9, c.grey, 'filled', ...
        'MarkerFaceAlpha',0.5, 'HandleVisibility','off');
end
scatter(NaN,NaN,9,c.grey,'filled','DisplayName','fixedPoints (independent check)');

xlabel('Load L (cmH_2O)'); ylabel('Equilibrium capacity C^* (cmH_2O)');
title(sprintf('Bistability at C_{max} = %g cmH_2O', C_max));
legend('Location','southwest','Box','off');
ylim(yl); xlim([0 br.L_fold_high*1.5]);
grid on; set(gca,'GridAlpha',0.08);

text(mean([br.L_fold_low br.L_fold_high]), C_max*0.99, 'bistable', ...
    'HorizontalAlignment','center','FontSize',9,'Color',c.fold);

% ---------------- Panel 2: hysteresis ----------------
nexttile; hold on;

Lup   = linspace(0, br.L_fold_high*1.3, 400);
Ldown = flip(Lup);

Cup   = sweepLoad(Lup,   C_max, p, C_max*0.999);
Cdown = sweepLoad(Ldown, C_max, p, C_max*p.alpha/(p.alpha+p.beta)*1.001);

plot(Lup,   Cup,   '-',  'Color', c.ETT,   'DisplayName','L increasing (loading)');
plot(Ldown, Cdown, '-',  'Color', c.TRACH, 'DisplayName','L decreasing (unloading)');

xline(br.L_fold_low,  ':', 'Color',c.fold, 'HandleVisibility','off');
xline(br.L_fold_high, ':', 'Color',c.fold, 'HandleVisibility','off');

xlabel('Load L (cmH_2O)'); ylabel('Capacity C (cmH_2O)');
title('Hysteresis: the path back is not the path out');
legend('Location','southwest','Box','off');
ylim(yl); xlim([0 br.L_fold_high*1.5]);
grid on; set(gca,'GridAlpha',0.08);

text(br.L_fold_high, C_max*0.5, ' collapse', 'Color',c.fold,'FontSize',9);
text(br.L_fold_low,  C_max*0.2, ' recovery ', 'Color',c.fold,'FontSize',9, ...
    'HorizontalAlignment','right');

title(tl, 'H1: the load-capacity system is bistable and hysteretic', ...
    'FontWeight','bold','FontSize',13);

viz.save(fig, 'fig_bifurcation_1p');

% ---------------- Table ----------------
nf = lc.normalizedFolds(p);
T = table({'left';'right'}, [br.L_fold_low; br.L_fold_high], ...
          [nf.l_low; nf.l_high], [nf.x_low*C_max; nf.x_high*C_max], ...
    'VariableNames', {'fold','L_cmH2O','l_normalised','C_at_fold'});
writetable(T, fullfile(wob.projectRoot,'results','tables','m2_exp1_folds.csv'));

fprintf('\n  Folds at C_max=%g: L = [%.3f, %.3f] cmH2O; normalised l = [%.4f, %.4f]\n', ...
    C_max, br.L_fold_low, br.L_fold_high, nf.l_low, nf.l_high);
fprintf('  Hysteresis width: %.3f cmH2O of load\n', br.L_fold_high - br.L_fold_low);

% Sanity: inside the window there must be 3 fixed points; outside, 1.
Lin  = mean([br.L_fold_low br.L_fold_high]);
Lout = br.L_fold_high * 1.2;
fprintf('  Inside window  (L=%.2f): %d fixed points %s\n', Lin, ...
    numel(lc.fixedPoints(Lin, C_max, p).C), ternary(numel(lc.fixedPoints(Lin,C_max,p).C)==3,'PASS','FAIL'));
fprintf('  Outside window (L=%.2f): %d fixed points %s\n', Lout, ...
    numel(lc.fixedPoints(Lout, C_max, p).C), ternary(numel(lc.fixedPoints(Lout,C_max,p).C)==1,'PASS','FAIL'));
end

% ---------------- helpers ----------------

function C = sweepLoad(Lseq, C_max, p, C0)
% Quasi-static sweep: at each load, relax the state a long way towards its
% attractor. Long enough that the trajectory tracks the branch except where
% the branch ceases to exist -- which is exactly where the jump happens.
C = zeros(size(Lseq));
state = C0;
for k = 1:numel(Lseq)
    [~, y] = ode45(@(t,CC) lc.dCdt(max(CC,1e-6), Lseq(k), C_max, p), [0 40], state);
    state = y(end);
    C(k) = state;
end
end

function plotRuns(x, y, mask, col, ls, lw, name)
% Plot only the samples in `mask`, breaking the line where the mask does, so
% a dashed unstable arm never joins across a gap.
d = diff([false, mask(:)', false]);
starts = find(d == 1);
stops  = find(d == -1) - 1;
first = true;
for k = 1:numel(starts)
    idx = starts(k):stops(k);
    plot(x(idx), y(idx), ls, 'Color', col, 'LineWidth', lw, ...
        'HandleVisibility', ternary(first,'on','off'), 'DisplayName', name);
    first = false;
end
end

function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end
