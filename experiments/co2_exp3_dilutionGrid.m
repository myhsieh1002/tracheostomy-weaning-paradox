function T = co2_exp3_dilutionGrid(cfgCO2, cfgM1)
%CO2_EXP3_DILUTIONGRID H3 -- the dead-space benefit is diluted. THESIS FIGURE.
%
%   T = co2_exp3_dilutionGrid() maps the share of total ventilatory demand
%   that the dead-space saving accounts for, across (VCO2 x V_D/V_T).
%
%   H3, AS REDEFINED
%   ----------------
%   The build spec's H3 reads as "the dead-space benefit's share of total
%   ventilatory demand, diluted by VCO2 and V_D/V_T". Its section 2.2 wrote
%   the benefit as a device term, dV_D_bypassed(device) -- but Model 1
%   established that the bypass CANCELS between an ETT and a tracheostomy, so
%   read that way H3 would just restate Model 1's finding on a new axis label.
%
%   So both contrasts are reported, because they answer different questions:
%
%     * ETT -> trach       the apparatus difference. Small, and it is what
%                          a clinician actually chooses between.
%     * not intubated ->   the upper-airway bypass. Large, and it is what
%       trach              the 72 mL literature is about (Nunn, Chadda).
%
%   The DILUTION is the point, and it is the gas-exchange analogue of Model
%   1's f_device dilution: as CO2 production rises and V/Q mismatch worsens,
%   total ventilatory demand grows while the dead-space saving does not, so
%   its share shrinks. Two independent axes -- mechanics and gas exchange --
%   both say the device moves less exactly where disease is worse.

arguments
    cfgCO2 (1,1) struct = co2.loadConfig()
    cfgM1  (1,1) struct = wob.loadConfig()
end

c = viz.style();

VCO2s  = cfgCO2.metabolism.VCO2_grid;
VDVTs  = cfgCO2.deadspace.VD_VT_grid;
armN   = cfgCO2.arms.native;
armE   = cfgCO2.arms.ett;
armT   = cfgCO2.arms.trach;

nV = numel(VCO2s); nD = numel(VDVTs);
fracDevice = nan(nD, nV);   % ETT -> trach
fracBypass = nan(nD, nV);   % not intubated -> trach
rows = {};

for i = 1:nD
    cfgD = cfgCO2;
    cfgD.deadspace.VD_VT = VDVTs(i);
    phi = co2.alveolarDeadSpaceFraction(cfgD, cfgM1);

    for j = 1:nV
        vN = co2.requiredVentilation(cfgD, cfgM1, armN, VCO2=VCO2s(j), Phi=phi);
        vE = co2.requiredVentilation(cfgD, cfgM1, armE, VCO2=VCO2s(j), Phi=phi);
        vT = co2.requiredVentilation(cfgD, cfgM1, armT, VCO2=VCO2s(j), Phi=phi);

        dDevice = vE.V_E - vT.V_E;
        dBypass = vN.V_E - vT.V_E;

        fracDevice(i,j) = dDevice / vE.V_E;
        fracBypass(i,j) = dBypass / vN.V_E;

        rows(end+1,:) = {VCO2s(j), VDVTs(i), phi, vN.V_E, vE.V_E, vT.V_E, ...
            dDevice, dBypass, fracDevice(i,j), fracBypass(i,j), ...
            vE.V_T, vT.V_T, vE.ds.total*1e3, vT.ds.total*1e3}; %#ok<AGROW>
    end
end

T = cell2table(rows, 'VariableNames', ...
    {'VCO2','VD_VT','phi','VE_native','VE_ETT','VE_trach','dVE_device','dVE_bypass', ...
     'frac_device','frac_bypass','V_T_ETT','V_T_trach','Vd_ETT_mL','Vd_trach_mL'});
writetable(T, fullfile(wob.projectRoot,'results','tables','co2_exp3_dilution.csv'));

% ================= Figure =================
fig = figure('Position',[60 60 1280 440]);
tl = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

maps = {fracBypass, fracDevice};
titles = {'not intubated \rightarrow trach', 'ETT \rightarrow trach'};
subs   = {'the upper-airway bypass (what the 72 mL literature describes)', ...
          'the apparatus difference (what a clinician chooses between)'};

for m = 1:2
    nexttile; hold on;
    M = maps{m};
    imagesc(VCO2s, VDVTs, M); set(gca,'YDir','normal');
    colormap(gca, viridisLike());
    axis tight;
    xticks(VCO2s); yticks(VDVTs);
    xlabel('V̇CO_2 (mL/min)  \rightarrow  more catabolic');
    if m == 1, ylabel('V_D/V_T  \rightarrow  worse V/Q mismatch'); end
    title(titles{m});
    subtitle(subs{m}, 'FontSize',8, 'Color', c.axisGrey);
    cb = colorbar; cb.Label.String = 'share of ventilatory demand';
    for i = 1:nD
        for j = 1:nV
            text(VCO2s(j), VDVTs(i), sprintf('%.3f', M(i,j)), ...
                'HorizontalAlignment','center','FontSize',8, ...
                'Color', ternary(M(i,j) > 0.5*max(M(:)), [0 0 0], [1 1 1]));
        end
    end
end

% --- dilution curves ---
nexttile; hold on;
for i = 1:nD
    shade = 0.75 * (nD - i) / max(nD-1,1);
    colB = c.TRACH + (1 - c.TRACH)*shade;
    colD = c.ETT   + (1 - c.ETT)*shade;
    plot(VCO2s, fracBypass(i,:), '-o', 'Color', colB, 'MarkerFaceColor', colB, ...
        'MarkerSize',3.5, 'HandleVisibility', ternary(i==nD,'on','off'), ...
        'DisplayName','not intubated \rightarrow trach');
    plot(VCO2s, fracDevice(i,:), '-s', 'Color', colD, 'MarkerFaceColor', colD, ...
        'MarkerSize',3.5, 'HandleVisibility', ternary(i==nD,'on','off'), ...
        'DisplayName','ETT \rightarrow trach');
end
xlabel('V̇CO_2 (mL/min)'); ylabel('share of ventilatory demand');
title({'Dilution on both contrasts','(darker = worse V/Q mismatch)'});
legend('Location','northeast','Box','off');
grid on; set(gca,'GridAlpha',0.08);
set(gca,'YScale','log');

title(tl, 'H3: the dead-space saving is diluted by metabolic rate and V/Q mismatch — the gas-exchange analogue of Model 1''s f_{device}', ...
    'FontWeight','bold','FontSize',12);

viz.save(fig, 'fig_co2_dilution_grid');

% ================= Verdict =================
fprintf('\n  Dilution check (share must fall as VCO2 rises and as V_D/V_T rises):\n');
for nm = ["frac_bypass","frac_device"]
    M = ternary(nm=="frac_bypass", fracBypass, fracDevice);
    monoV = all(diff(M,1,2) < 0, 'all');
    monoD = all(diff(M,1,1) < 0, 'all');
    fprintf('    %-12s falls with VCO2: %s | falls with V_D/V_T: %s | range %.4f - %.4f\n', ...
        nm, ternary(monoV,'PASS','FAIL'), ternary(monoD,'PASS','FAIL'), min(M(:)), max(M(:)));
end

fprintf('\n  The two contrasts, at the base case (VCO2=%g, V_D/V_T=%.2f):\n', VCO2s(1), VDVTs(1));
sel = find(T.VCO2==VCO2s(1) & T.VD_VT==VDVTs(1), 1);
fprintf('    not intubated -> trach: dV_E = %.3f L/min (%.1f%% of demand)\n', ...
    T.dVE_bypass(sel), 100*T.frac_bypass(sel));
fprintf('    ETT -> trach:           dV_E = %.3f L/min (%.1f%% of demand)  <- the device step\n', ...
    T.dVE_device(sel), 100*T.frac_device(sel));
fprintf('    ratio: the device step is %.0f%% of the bypass step\n', ...
    100*T.dVE_device(sel)/T.dVE_bypass(sel));

fprintf('\n  Worst case (VCO2=%g, V_D/V_T=%.2f):\n', VCO2s(end), VDVTs(end));
sel2 = find(T.VCO2==VCO2s(end) & T.VD_VT==VDVTs(end), 1);
fprintf('    ETT -> trach saves %.3f L/min = %.1f%% of a %.1f L/min demand\n', ...
    T.dVE_device(sel2), 100*T.frac_device(sel2), T.VE_ETT(sel2));
end

function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end

function m = viridisLike()
anchors = [0.267 0.005 0.329; 0.283 0.141 0.458; 0.254 0.265 0.530;
           0.207 0.372 0.553; 0.164 0.471 0.558; 0.128 0.567 0.551;
           0.135 0.659 0.518; 0.267 0.749 0.441; 0.478 0.821 0.318;
           0.741 0.873 0.150; 0.993 0.906 0.144];
x = linspace(0,1,size(anchors,1)); xi = linspace(0,1,256);
m = [interp1(x,anchors(:,1),xi)', interp1(x,anchors(:,2),xi)', interp1(x,anchors(:,3),xi)'];
end
