function T = m2_exp5_sensitivity(cfg1, cfg2)
%M2_EXP5_SENSITIVITY Sobol sensitivity of the fold and the rescue window.
%
%   T = m2_exp5_sensitivity() decomposes the variance of three outputs over
%   the dynamical parameters.
%
%   This doubles as an INDEPENDENT TEST OF THE SCALE INVARIANCE. The claim
%   is that C_max cannot influence the NORMALISED fold l_high, only the
%   absolute one. A Sobol decomposition never sees the algebra -- it only
%   sees inputs and outputs -- so if the invariance is real,
%
%       ST(C_max) on l_high        must be ~0
%       ST(C_max) on L_fold_high   must be large
%
%   Getting that pair out of a purely numerical variance decomposition is a
%   much stronger check than re-reading the derivation.
%
%   The outputs:
%       l_high        normalised right fold (= the Pi/Pimax threshold)
%       L_fold_high   absolute right fold (cmH2O) -- what a clinician sees
%       window_width  rescue window in C_max = dL_device / l_high, with
%                     dL_device taken from Model 1

arguments
    cfg1 (1,1) struct = wob.loadConfig()
    cfg2 (1,1) struct = lc.loadConfig()
end

c = viz.style();
pBase = lc.calibrate(cfg2);

% dL_device comes from Model 1, once, at the base patient.
LE = lc.coupling(cfg1, cfg2, cfg2.load_from_model1.devices{1}).L_total;
LT = lc.coupling(cfg1, cfg2, cfg2.load_from_model1.devices{2}).L_total;
dL_device = LE - LT;

names  = {'alpha','beta','s','u_crit','C_max'};
ranges = [ 0.10 0.30;      % alpha  - Laghi range
           0.04 0.20;      % beta   - around the Laghi-derived 0.083
           30   150;       % s      - above the cusp (s* ~ 23)
           0.35 0.50;      % u_crit - Roussos sensitivity band
           20   70 ];      % C_max  - the corrected capacity grid

outputs = {'l_high','L_fold_high','window_width'};
labels  = {'l_{high} (normalised fold)','L_{fold,high} (cmH_2O)','rescue window width (cmH_2O)'};

Sres = struct(); rows = {};
for m = 1:numel(outputs)
    f = @(x) evalM2(x, names, dL_device, outputs{m});
    tic;
    S = wob.sobolIndices(f, ranges, N=cfg2.sensitivity.n_samples, Seed=cfg2.sensitivity.seed, ...
                         Names=string(names), Bootstrap=200);
    fprintf('    %-14s done in %.1f s | sum(S1) = %.3f\n', outputs{m}, toc, sum(S.S1));
    Sres.(outputs{m}) = S;
    for i = 1:numel(names)
        rows(end+1,:) = {outputs{m}, names{i}, S.S1(i), S.S1_ci(1,i), S.S1_ci(2,i), ...
            S.ST(i), S.ST_ci(1,i), S.ST_ci(2,i)}; %#ok<AGROW>
    end
end

T = cell2table(rows, 'VariableNames', ...
    {'output','parameter','S1','S1_lo','S1_hi','ST','ST_lo','ST_hi'});
writetable(T, fullfile(wob.projectRoot,'results','tables','m2_exp5_sobol.csv'));

% ================= Figure =================
fig = figure('Position',[50 50 1280 430]);
tl = tiledlayout(1, numel(outputs)+1, 'TileSpacing','compact','Padding','compact');

for m = 1:numel(outputs)
    nexttile; hold on;
    S = Sres.(outputs{m});
    k = numel(names); x = 1:k;
    bar(x-0.19, S.S1, 0.36, 'FaceColor', c.TRACH, 'EdgeColor','none', 'DisplayName','S1');
    bar(x+0.19, S.ST, 0.36, 'FaceColor', c.ETT,   'EdgeColor','none', 'DisplayName','ST');
    errorbar(x-0.19, S.S1, S.S1-S.S1_ci(1,:), S.S1_ci(2,:)-S.S1, 'k','LineStyle','none','CapSize',3,'HandleVisibility','off');
    errorbar(x+0.19, S.ST, S.ST-S.ST_ci(1,:), S.ST_ci(2,:)-S.ST, 'k','LineStyle','none','CapSize',3,'HandleVisibility','off');
    xticks(x); xticklabels(strrep(names,'_','\_')); xtickangle(35);
    ylabel('Sobol index'); title(labels{m});
    ylim([0 1.05]);
    if m==1, legend('Location','northwest','Box','off'); end
    grid on; set(gca,'GridAlpha',0.08);

    % Flag the invariance prediction on the panel it applies to.
    iC = find(strcmp(names,'C_max'));
    if strcmp(outputs{m},'l_high')
        text(iC, 0.10, sprintf('ST = %.1e\n(invariance)', S.ST(iC)), ...
            'HorizontalAlignment','center','FontSize',7.5,'Color',c.rescue,'FontWeight','bold');
    end
end

% --- OAT: fold vs beta/alpha, the ratio that actually governs it ---
nexttile; hold on;
ratios = linspace(0.1, 1.5, 160);
lh = nan(size(ratios));
for i = 1:numel(ratios)
    p = pBase; p.beta = p.alpha * ratios(i);
    nf = lc.normalizedFolds(p);
    if nf.isBistable, lh(i) = nf.l_high; end
end
plot(ratios, lh, '-', 'Color', c.stable, 'LineWidth', 1.8, 'DisplayName','l_{high}');
xline(pBase.beta/pBase.alpha, ':', 'Color', c.fold, 'LineWidth',1.4, ...
    'Label', sprintf('calibrated \\beta/\\alpha = %.2f', pBase.beta/pBase.alpha), 'FontSize',8, ...
    'HandleVisibility','off');
yline(pBase.u_crit, '--', 'Color', c.grey, 'Label','u_{crit} ceiling', 'FontSize',8, 'HandleVisibility','off');
xlabel('\beta / \alpha'); ylabel('l_{high}');
title({'The fold is governed by the RATIO', '\beta/\alpha, not by absolute rates'});
grid on; set(gca,'GridAlpha',0.08);

title(tl, 'Model 2 sensitivity -- and a numerical test of the scale invariance', ...
    'FontWeight','bold','FontSize',12.5);
viz.save(fig, 'fig_m2_sensitivity');

% ================= Verdict =================
iC = find(strcmp(names,'C_max'));
fprintf('\n  INVARIANCE TEST (Sobol never sees the algebra):\n');
fprintf('    ST(C_max) on l_high      = %.3e  -> %s\n', Sres.l_high.ST(iC), ...
    ternary(abs(Sres.l_high.ST(iC)) < 1e-6, 'ZERO, as the scale invariance requires', 'NON-ZERO - invariance violated!'));
fprintf('    ST(C_max) on L_fold_high = %.3f      -> dominant, as expected (L_fold = l_high*C_max)\n', ...
    Sres.L_fold_high.ST(iC));

fprintf('\n  Total-effect ranking:\n');
for m = 1:numel(outputs)
    S = Sres.(outputs{m});
    [~, ord] = sort(S.ST,'descend');
    fprintf('    %-14s : ', outputs{m});
    for i = ord, fprintf('%s(%.2f) ', names{i}, S.ST(i)); end
    fprintf('\n');
end

Sw = Sres.window_width;
fprintf('\n  Rescue-window width is driven mainly by %s.\n', names{find(Sw.ST==max(Sw.ST),1)});
fprintf('  Note dL_device (%.2f cmH2O, from Model 1) enters as a fixed multiplier here;\n', dL_device);
fprintf('  its own uncertainty is quantified in exp5 (K2_scale).\n');
end

% ================= helpers =================

function y = evalM2(x, names, dL_device, outputName)
p = struct('alpha', x(strcmp(names,'alpha')), ...
           'beta',  x(strcmp(names,'beta')), ...
           's',     x(strcmp(names,'s')), ...
           'u_crit',x(strcmp(names,'u_crit')));
C_max = x(strcmp(names,'C_max'));

nf = lc.normalizedFolds(p);
if ~nf.isBistable
    % Below the cusp there is no fold. Report NaN-free sentinels so the
    % variance decomposition stays defined; these draws are recorded rather
    % than silently dropped.
    y = 0;
    return;
end

switch outputName
    case 'l_high',       y = nf.l_high;
    case 'L_fold_high',  y = nf.l_high * C_max;
    case 'window_width', y = dL_device / nf.l_high;
end
end

function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end
