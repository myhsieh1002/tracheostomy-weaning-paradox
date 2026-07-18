function T = m2_exp4_earlyWarning(cfg2)
%M2_EXP4_EARLYWARNING H4 -- critical slowing down as a weaning early-warning signal.
%
%   T = m2_exp4_earlyWarning() ramps the load slowly towards the right fold
%   under noise and tracks the rolling variance and lag-1 autocorrelation of
%   the capacity.
%
%   H4: approaching a saddle-node, the restoring eigenvalue goes to zero, so
%   perturbations decay ever more slowly. Variance and AR(1) rise. Clinically
%   this maps onto rising breathing variability before weaning failure.
%
%   The test is a monotone-trend statistic (Kendall's tau) against distance
%   to the fold, over an ENSEMBLE of noise realisations -- a single path is
%   not evidence, since these statistics are themselves noisy.
%
%   THE HONEST CAVEAT
%   -----------------
%   This is where alpha's absolute value finally matters. Everything in
%   H1-H3 depends only on beta/alpha, but the RATE of critical slowing down
%   is set by the absolute eigenvalue, and alpha is only semi-constrained
%   (Laghi 1995's recovery is biphasic and no single alpha fits it). So the
%   qualitative rise is a robust prediction of the fold; the TIMESCALE over
%   which it becomes detectable is not, and should not be quoted as a lead
%   time.

arguments
    cfg2 (1,1) struct = lc.loadConfig()
end

c = viz.style();
p = lc.calibrate(cfg2);
nf = lc.normalizedFolds(p);

C_max  = 50;
L_fold = nf.l_high * C_max;

% Ramp up to just short of the fold, slowly enough that the state tracks the
% branch quasi-statically -- otherwise we would be watching the ramp outrun
% the dynamics rather than critical slowing down.
T_end  = 1200;                     % hours
nSteps = 240000;                   % dt = 0.005 h, for integration accuracy
t   = linspace(0, T_end, nSteps);
L0  = nf.l_low * C_max * 1.02;
L1  = L_fold * 0.99;
L_t = L0 + (L1 - L0) * (t / T_end);

% The statistics are computed on a SUBSAMPLED series. The integrator needs a
% small step; AR(1) needs a sampling interval comparable to the relaxation
% time 1/alpha (~6.7 h here), or it saturates at ~1 and reports nothing.
subN  = round(2.0 / (T_end/nSteps));   % -> ~2 h between analysed samples
winN  = 100;                           % -> ~200 h window

nEns  = 24;
sigma = cfg2.dynamics.noise_sigma;

% The separatrix at the START of the ramp: any path below it has escaped to
% the low branch. Statistics computed after that are describing a collapsed
% system, not an approach to the fold, so each path is truncated there.
fp0 = lc.fixedPoints(L0, C_max, p);
separatrix = fp0.C(2);

varAll = []; ar1All = []; tMid = [];
paths = zeros(nEns, nSteps);
collapseT = nan(1, nEns);
for e = 1:nEns
    s = lc.simulateStochastic(t, L_t, C_max, p, Sigma=sigma, Seed=100+e);
    paths(e,:) = s.C;

    iCollapse = find(s.C < separatrix, 1, 'first');
    if ~isempty(iCollapse)
        collapseT(e) = t(iCollapse);
    end

    ew = lc.earlyWarning(t, s.C, Window=winN, Subsample=subN);
    if isempty(varAll)
        tMid = ew.tMid;
        varAll = nan(nEns, numel(tMid));
        ar1All = nan(nEns, numel(tMid));
    end
    v = ew.variance; a = ew.ar1;
    if ~isnan(collapseT(e))
        v(tMid >= collapseT(e)) = NaN;   % drop post-collapse windows
        a(tMid >= collapseT(e)) = NaN;
    end
    varAll(e,:) = v; %#ok<AGROW>
    ar1All(e,:) = a; %#ok<AGROW>
end

% Keep only windows where most of the ensemble is still pre-collapse;
% averaging over a handful of survivors is a survivorship artefact.
enough = sum(~isnan(varAll), 1) >= 0.6*nEns;
tMid = tMid(enough); varAll = varAll(:,enough); ar1All = ar1All(:,enough);

% Distance to the fold at each window centre.
L_mid = interp1(t, L_t, tMid);
dist  = L_fold - L_mid;

varMean = mean(varAll,1,'omitnan'); varLo = prctile(varAll,10,1); varHi = prctile(varAll,90,1);
ar1Mean = mean(ar1All,1,'omitnan'); ar1Lo = prctile(ar1All,10,1); ar1Hi = prctile(ar1All,90,1);

% Monotone trend vs distance-to-fold. Kendall's tau, not Pearson: we are
% testing for a monotone rise, not a linear one, and these statistics are
% heavy-tailed.
tauVar = corr(dist', varMean', 'Type','Kendall');
tauAr1 = corr(dist', ar1Mean', 'Type','Kendall');

% ================= Figure =================
fig = figure('Position',[60 60 1240 700]);
tl = tiledlayout(3,2,'TileSpacing','compact','Padding','compact');

% --- capacity paths ---
nexttile([1 2]); hold on;
for e = 1:min(nEns,8)
    plot(t, paths(e,:), '-', 'Color',[c.grey 0.5], 'LineWidth',0.5, 'HandleVisibility','off');
end
plot(t, mean(paths,1), '-', 'Color',c.stable, 'LineWidth',1.6, 'DisplayName','ensemble mean');
yyaxis right;
plot(t, L_t, '--', 'Color',c.ETT, 'LineWidth',1.4, 'DisplayName','load L(t)');
yline(L_fold, ':', 'Color',c.fold, 'LineWidth',1.4, 'Label','right fold', 'FontSize',8, ...
    'HandleVisibility','off');
ylabel('Load L (cmH_2O)'); set(gca,'YColor',c.ETT);
yyaxis left; ylabel('Capacity C (cmH_2O)'); set(gca,'YColor',c.stable);
xlabel('Time (hours)');
title(sprintf('Slow ramp towards the fold (%d noise realisations, \\sigma = %.3f\\cdotC_{max})', nEns, sigma));
legend('Location','southwest','Box','off');
grid on; set(gca,'GridAlpha',0.08);

% --- variance vs time ---
nexttile; hold on;
fill([tMid fliplr(tMid)], [varLo fliplr(varHi)], c.TRACH, 'FaceAlpha',0.2, 'EdgeColor','none');
plot(tMid, varMean, '-', 'Color',c.TRACH, 'LineWidth',1.6);
xlabel('Time (hours)'); ylabel('Rolling variance (detrended)');
title('Variance rises as the fold nears');
grid on; set(gca,'GridAlpha',0.08);

% --- AR(1) vs time ---
nexttile; hold on;
fill([tMid fliplr(tMid)], [ar1Lo fliplr(ar1Hi)], c.ETT, 'FaceAlpha',0.2, 'EdgeColor','none');
plot(tMid, ar1Mean, '-', 'Color',c.ETT, 'LineWidth',1.6);
xlabel('Time (hours)'); ylabel('Lag-1 autocorrelation');
title('AR(1) rises as the fold nears');
grid on; set(gca,'GridAlpha',0.08);

% --- vs distance to fold ---
nexttile; hold on;
plot(dist, varMean, '-', 'Color',c.TRACH, 'LineWidth',1.6);
set(gca,'XDir','reverse');
xlabel('Distance to fold, L_{fold} - L (cmH_2O)   \rightarrow approaching');
ylabel('Rolling variance');
title(sprintf('Kendall \\tau = %.2f', tauVar));
grid on; set(gca,'GridAlpha',0.08);

nexttile; hold on;
plot(dist, ar1Mean, '-', 'Color',c.ETT, 'LineWidth',1.6);
set(gca,'XDir','reverse');
xlabel('Distance to fold, L_{fold} - L (cmH_2O)   \rightarrow approaching');
ylabel('Lag-1 autocorrelation');
title(sprintf('Kendall \\tau = %.2f', tauAr1));
grid on; set(gca,'GridAlpha',0.08);

title(tl, 'H4: critical slowing down before the fold -- a mechanistic basis for rising breathing variability', ...
    'FontWeight','bold','FontSize',12.5);

viz.save(fig, 'fig_early_warning');

% ================= Table =================
T = table(tMid', L_mid', dist', varMean', ar1Mean', ...
    'VariableNames', {'t_hr','L','dist_to_fold','variance','ar1'});
writetable(T, fullfile(wob.projectRoot,'results','tables','m2_exp4_early_warning.csv'));

% ================= Verdict =================
fprintf('\n  Ensemble: %d realisations, sigma = %.3f*C_max, C_max = %g\n', nEns, sigma, C_max);
fprintf('  Fold at L = %.3f cmH2O; ramp %.2f -> %.2f over %g h\n', L_fold, L0, L1, T_end);
fprintf('  Analysed sampling: every %.1f h (relaxation time 1/alpha = %.1f h); window %.0f h\n', ...
    subN*(T_end/nSteps), 1/p.alpha, winN*subN*(T_end/nSteps));
fprintf('  Noise-induced collapse before the fold in %d/%d paths (median t = %.0f h)\n', ...
    sum(~isnan(collapseT)), nEns, median(collapseT,'omitnan'));
fprintf('  Post-collapse windows excluded; analysis kept where >=60%% of paths survive.\n');
fprintf('  Variance: %.3g (far) -> %.3g (near fold), x%.1f\n', ...
    varMean(1), varMean(end), varMean(end)/varMean(1));
fprintf('  AR(1)   : %.3f (far) -> %.3f (near fold)\n', ar1Mean(1), ar1Mean(end));
fprintf('  Monotone trend vs distance-to-fold (Kendall tau, negative = rises on approach):\n');
fprintf('    variance tau = %+.3f -> %s\n', tauVar, ternary(tauVar < -0.5,'PASS','WEAK/FAIL'));
fprintf('    AR(1)    tau = %+.3f -> %s\n', tauAr1, ternary(tauAr1 < -0.5,'PASS','WEAK/FAIL'));
fprintf('\n  Caveat: the RATE of slowing is set by the absolute eigenvalue, hence by alpha,\n');
fprintf('  which is only semi-constrained. The qualitative rise is robust; the lead time is not.\n');
end

function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end
