function T = m2_exp2_cusp2p(cfg2)
%M2_EXP2_CUSP2P H2 -- two-parameter structure. The cusp is NOT where the spec put it.
%
%   T = m2_exp2_cusp2p() traces the saddle-node loci in two parameter
%   planes and locates the cusp point.
%
%   WHAT THE SPEC EXPECTED, AND WHAT THE MODEL DOES
%   -----------------------------------------------
%   The spec predicted a cusp in the (L, C_max) plane. There is none, and
%   the reason is the scale invariance: the folds sit at fixed NORMALISED
%   load, so in absolute load they lie on
%
%       L_fold = l_fold * C_max
%
%   -- exactly straight lines through the origin (verified to a residual of
%   ~1e-14). Two straight lines through the origin meet only at the origin.
%   The bistable region is therefore a WEDGE with its apex at C_max = 0, not
%   a cusp. The spec's three regions are all present; only the geometry, and
%   what it implies, is different.
%
%   This matters for the argument rather than being a technicality. A cusp
%   would mean bistability disappears at some finite capacity -- that below
%   a critical C_max the tipping point simply ceases to exist. It does not.
%   The bistable window persists at every capacity, shrinking in absolute
%   load but never closing. What changes with C_max is only WHERE the fold
%   sits, which is why the rescue window moves rather than vanishes (see
%   m2_exp3).
%
%   WHERE THE CUSP ACTUALLY IS
%   --------------------------
%   Bistability has to be born somewhere. Since the fold structure depends
%   only on (beta/alpha, s, u_crit), the cusp lives in those coordinates,
%   not in C_max. Panel B sweeps s: below a critical steepness the two folds
%   merge and the system is monostable at every load. THAT is the cusp.

arguments
    cfg2 (1,1) struct = lc.loadConfig()
end

c = viz.style();
p = lc.calibrate(cfg2);

fig = figure('Position',[60 60 1200 460]);
tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

% ================= Panel A: (L, C_max) -- the wedge =================
nexttile; hold on;

nf = lc.normalizedFolds(p);
CmaxAxis = linspace(0, 80, 200);
L_low  = nf.l_low  * CmaxAxis;
L_high = nf.l_high * CmaxAxis;

Lmax = 32;
% monostable-sustainable (below the left fold)
patch([0 L_low  fliplr([0 0])], [0 CmaxAxis 80 0], c.TRACH_light, ...
    'EdgeColor','none', 'FaceAlpha',0.75, 'DisplayName','monostable: sustainable');
% bistable wedge
patch([L_low fliplr(L_high)], [CmaxAxis fliplr(CmaxAxis)], c.rescue, ...
    'EdgeColor','none', 'FaceAlpha',0.55, 'DisplayName','bistable');
% monostable-failure (above the right fold)
patch([L_high fliplr([Lmax Lmax])], [CmaxAxis 80 0], c.ETT_light, ...
    'EdgeColor','none', 'FaceAlpha',0.75, 'DisplayName','monostable: failure');

plot(L_low,  CmaxAxis, '-', 'Color', c.fold, 'LineWidth', 1.8, 'DisplayName','saddle-node loci');
plot(L_high, CmaxAxis, '-', 'Color', c.fold, 'LineWidth', 1.8, 'HandleVisibility','off');
plot(0, 0, 'o', 'MarkerSize',7, 'MarkerFaceColor',c.fold, 'MarkerEdgeColor','w', ...
    'DisplayName','apex (C_{max} = 0)');

xlabel('Load L (cmH_2O)'); ylabel('C_{max} (cmH_2O)');
title('(L, C_{max}): a WEDGE, not a cusp');
subtitle(sprintf('loci are exactly L = %.4f\\cdotC_{max} and %.4f\\cdotC_{max}', nf.l_low, nf.l_high), ...
    'FontSize',9, 'Color', c.axisGrey);
legend('Location','southeast','Box','off');
xlim([0 Lmax]); ylim([0 80]);
grid on; set(gca,'GridAlpha',0.08);

% ================= Panel B: (l, s) -- the real cusp =================
nexttile; hold on;

sAxis = linspace(15, 120, 400);
lLow = nan(size(sAxis)); lHigh = nan(size(sAxis));
for k = 1:numel(sAxis)
    pk = p; pk.s = sAxis(k);
    nfk = lc.normalizedFolds(pk);
    if nfk.isBistable
        lLow(k) = nfk.l_low; lHigh(k) = nfk.l_high;
    end
end

sCusp = findCusp(p);
pC = p; pC.s = sCusp;
nfC = lc.normalizedFolds(pC);
lCusp = mean([nfC.l_low, nfC.l_high]);

ok = ~isnan(lLow);
patch([lLow(ok) fliplr(lHigh(ok))], [sAxis(ok) fliplr(sAxis(ok))], c.rescue, ...
    'EdgeColor','none','FaceAlpha',0.5, 'DisplayName','bistable');
plot(lLow(ok),  sAxis(ok), '-', 'Color',c.fold, 'LineWidth',1.8, 'DisplayName','saddle-node loci');
plot(lHigh(ok), sAxis(ok), '-', 'Color',c.fold, 'LineWidth',1.8, 'HandleVisibility','off');

plot(lCusp, sCusp, 'p', 'MarkerSize',15, 'MarkerFaceColor',c.fold, 'MarkerEdgeColor','w', ...
    'LineWidth',1, 'DisplayName',sprintf('CUSP at s = %.1f', sCusp));

xline(p.u_crit, '--', 'Color', c.axisGrey, 'LineWidth',1.2, ...
    'Label','u_{crit} = 0.40 (ceiling as s\rightarrow\infty)', 'FontSize',8, ...
    'LabelVerticalAlignment','bottom', 'HandleVisibility','off');
yline(p.s, ':', 'Color', c.stable, 'LineWidth',1.2, ...
    'Label',sprintf('calibrated s = %.1f', p.s), 'FontSize',8, 'HandleVisibility','off');

xlabel('Normalised load l = L / C_{max}   ( = P_i/P_{imax} )');
ylabel('Sigmoid steepness s');
title('(l, s): here IS the cusp');
subtitle('below it the folds merge and bistability ceases to exist', ...
    'FontSize',9, 'Color', c.axisGrey);
legend('Location','southeast','Box','off');
grid on; set(gca,'GridAlpha',0.08);

title(tl, 'H2: three regions as predicted -- but the cusp lives in the fatigue-threshold sharpness, not in capacity', ...
    'FontWeight','bold','FontSize',12.5);

viz.save(fig, 'fig_cusp_2p');

% ================= Table =================
rows = {};
for k = 1:numel(sAxis)
    if ~isnan(lLow(k))
        rows(end+1,:) = {sAxis(k), lLow(k), lHigh(k), lHigh(k)-lLow(k)}; %#ok<AGROW>
    end
end
T = cell2table(rows, 'VariableNames', {'s','l_low','l_high','width'});
writetable(T, fullfile(wob.projectRoot,'results','tables','m2_exp2_cusp.csv'));

% ================= Reported numbers =================
fprintf('\n  (L, C_max) plane: loci are L = %.6f*C_max and %.6f*C_max\n', nf.l_low, nf.l_high);
resid = checkLinearity(p);
fprintf('    linearity residual over C_max in [5,160]: %.2e  -> straight through origin, NO cusp\n', resid);
fprintf('    => bistability NEVER disappears with falling capacity; the window MOVES.\n');
fprintf('\n  Cusp point (beta/alpha = %.3f, u_crit = %.2f): s* = %.2f, l* = %.4f\n', ...
    p.beta/p.alpha, p.u_crit, sCusp, lCusp);
fprintf('    calibrated s = %.1f is %s the cusp -> the model IS bistable\n', ...
    p.s, ternary(p.s > sCusp, 'ABOVE', 'BELOW'));
fprintf('    margin: s/s* = %.2f\n', p.s/sCusp);
end

% ================= helpers =================

function sCusp = findCusp(p)
%FINDCUSP Steepness at which the two folds merge.
%   Bisection on "is it bistable", which is a clean monotone predicate in s:
%   the fold width shrinks to zero from above as s falls.
lo = 1; hi = p.s;
if ~isBi(p, hi)
    error('m2_exp2:notBistable', 'The calibrated parameters are not bistable; nothing to bracket.');
end
for i = 1:80
    mid = 0.5*(lo+hi);
    if isBi(p, mid), hi = mid; else, lo = mid; end
end
sCusp = hi;
end

function tf = isBi(p, s)
p.s = s;
nf = lc.normalizedFolds(p);
tf = nf.isBistable;
end

function resid = checkLinearity(p)
Cm = [5 10 20 40 80 160];
Lh = zeros(size(Cm));
for k = 1:numel(Cm)
    br = lc.equilibriumBranch(Cm(k), p);
    Lh(k) = br.L_fold_high;
end
slope = Lh / Cm;                 % least-squares through the origin
resid = max(abs(Lh - slope*Cm));
end

function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end
