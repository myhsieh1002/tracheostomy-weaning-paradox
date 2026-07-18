function T = exp1_resistanceCurves(cfg)
%EXP1_RESISTANCECURVES H1 -- resistance vs flow for every device.
%
%   T = exp1_resistanceCurves() sweeps flow across each device, plots the
%   Rohrer pressure-drop curves, and checks the ID-scaling exponent against
%   the physics-based path.
%
%   H1: at matched ID, a tracheostomy tube has lower resistance than an ETT
%   (it is shorter), and both rise steeply as ID falls.

arguments
    cfg (1,1) struct = wob.loadConfig()
end

c = viz.style();
flow = linspace(cfg.experiments.exp1.flow_range_Lps(1), ...
                cfg.experiments.exp1.flow_range_Lps(2), ...
                cfg.experiments.exp1.n_flow);

% Real tubes only. NATIVE_UPPER_AIRWAY is a pseudo-device standing in for
% the intact airway (used by exp3) and has no internal diameter to plot.
names = fieldnames(cfg.devices);
names = names(~strcmp(names, 'NATIVE_UPPER_AIRWAY'));
rows = {};

fig = figure('Position', [100 100 1000 400]);

% --- Panel 1: pressure drop vs flow ---
subplot(1,2,1); hold on;
for k = 1:numel(names)
    dev = wob.getDevice(cfg, names{k});
    dP = wob.rohrerDrop(dev, flow);

    isETT = strcmp(dev.type, 'ETT');
    col = c.ETT * isETT + c.TRACH * ~isETT;
    shade = 1 - 0.35 * (dev.ID_mm - 7) / 1.0;
    plot(flow, dP, 'Color', min(col * shade, 1), ...
        'LineStyle', ternary(isETT, '-', '--'), 'DisplayName', strrep(names{k}, '_', ' '));

    rows(end+1,:) = {names{k}, dev.type, dev.ID_mm, dev.length_cm, dev.K1, dev.K2, ...
                     wob.rohrerDrop(dev, 1.0), wob.deviceResistance(dev, 1.0)}; %#ok<AGROW>
end
xlabel('Flow (L/s)'); ylabel('\DeltaP across device (cmH_2O)');
title('Device pressure drop (Rohrer)');
legend('Location','northwest','Box','off');
grid on; set(gca,'GridAlpha',0.08);

% --- Panel 2: ID scaling, empirical vs physics ---
subplot(1,2,2); hold on;
IDs = linspace(6, 9, 40);
refFlow = 1.0;

for L_cm = [27, 8]
    dPphys = arrayfun(@(d) abs(wob.physicsScaling(d, L_cm, refFlow).dP_total), IDs);
    isETTlen = L_cm > 15;
    plot(IDs, dPphys, 'Color', ternary(isETTlen, c.ETT, c.TRACH), ...
        'LineWidth', 1.2, 'DisplayName', sprintf('physics, L=%d cm', L_cm));

    % Local scaling exponent d(log dP)/d(log ID)
    expo = diff(log(dPphys)) ./ diff(log(IDs));
    fprintf('  physics L=%2d cm: dP ~ ID^%.2f at 1 L/s (expect -4 to -5)\n', L_cm, mean(expo));
end

for k = 1:numel(names)
    dev = wob.getDevice(cfg, names{k});
    isETT = strcmp(dev.type,'ETT');
    scatter(dev.ID_mm, wob.rohrerDrop(dev, refFlow), 45, ...
        ternary(isETT, c.ETT, c.TRACH), 'filled', ...
        'MarkerEdgeColor','w', 'HandleVisibility','off');
end
set(gca,'XScale','log','YScale','log');
xlabel('Internal diameter (mm)'); ylabel(sprintf('\\DeltaP at %.1f L/s (cmH_2O)', refFlow));
title('ID scaling: empirical points vs physics');
legend('Location','northeast','Box','off');
grid on; set(gca,'GridAlpha',0.08);

viz.save(fig, 'fig_resistance_curves');

T = cell2table(rows, 'VariableNames', ...
    {'device','type','ID_mm','length_cm','K1','K2','dP_at_1Lps','R_at_1Lps'});
writetable(T, fullfile(wob.projectRoot,'results','tables','exp1_resistance.csv'));

% --- H1 check: matched ID, trach < ETT ---
fprintf('\n  H1 check (matched ID 7.0 and 8.0, at 1 L/s):\n');
for ID = [7.0, 8.0]
    e = T(strcmp(T.type,'ETT')   & T.ID_mm==ID, :);
    t = T(strcmp(T.type,'TRACH') & T.ID_mm==ID, :);
    if ~isempty(e) && ~isempty(t)
        pass = t.dP_at_1Lps < e.dP_at_1Lps;
        fprintf('    ID %.1f: ETT %.2f vs TRACH %.2f cmH2O -> %s\n', ID, ...
            e.dP_at_1Lps, t.dP_at_1Lps, ternary(pass,'PASS','FAIL'));
    end
end
end

function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end
