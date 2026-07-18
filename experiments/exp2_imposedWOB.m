function T = exp2_imposedWOB(cfg)
%EXP2_IMPOSEDWOB H2 -- imposed work and PTP at matched ventilation (Mode A).
%
%   T = exp2_imposedWOB() compares every device under an IDENTICAL
%   prescribed waveform, so any difference is attributable to the tube
%   alone. Mode A asks nothing about the patient; it is the clean passive
%   contrast.
%
%   H2: at matched V_T and RR, the tracheostomy imposes less work and less
%   pressure-time product than the ETT.
%
%   Devices graded X (no published coefficients) are drawn hatched and
%   excluded from the H2 verdict.

arguments
    cfg (1,1) struct = wob.loadConfig()
end

c = viz.style();
devices = cellstr(cfg.experiments.exp2.devices);

rows = {};
for k = 1:numel(devices)
    dev = wob.getDevice(cfg, devices{k});
    r = wob.simulateModeA(cfg, dev);
    rows(end+1,:) = {devices{k}, dev.type, dev.ID_mm, dev.grade, ...
        r.WOB_device_J_L, r.WOB_device_J_min, r.PTP_device_per_min, ...
        r.peak_dP_device, r.mean_dP_device, r.peak_flow}; %#ok<AGROW>
end

T = cell2table(rows, 'VariableNames', ...
    {'device','type','ID_mm','grade','WOB_J_L','WOB_J_min','PTP_per_min', ...
     'peak_dP','mean_dP','peak_flow'});
writetable(T, fullfile(wob.projectRoot,'results','tables','exp2_table.csv'));

% ---------------- Figure ----------------
metrics = {'WOB_J_L','WOB_J_min','PTP_per_min'};
labels  = {'Imposed WOB (J/L)','Imposed WOB (J/min)','Imposed PTP (cmH_2O\cdots/min)'};

fig = figure('Position',[80 80 1160 400]);
tl = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

for m = 1:numel(metrics)
    nexttile; hold on;
    vals = T.(metrics{m});
    for k = 1:height(T)
        isETT = strcmp(T.type{k},'ETT');
        col = ternary(isETT, c.ETT, c.TRACH);
        unverified = strcmp(T.grade{k},'X');
        b = bar(k, vals(k), 0.68, 'FaceColor', col, 'EdgeColor','none');
        if unverified
            b.FaceAlpha = 0.35;
            b.EdgeColor = col; b.LineWidth = 1; b.LineStyle = ':';
        end
        text(k, vals(k), sprintf(' %.3g', vals(k)), 'HorizontalAlignment','center', ...
            'VerticalAlignment','bottom', 'FontSize',8);
    end
    xticks(1:height(T));
    xticklabels(strrep(T.device,'_',' '));
    xtickangle(35);
    ylabel(labels{m});
    ylim([0, max(vals)*1.18]);
    grid on; set(gca,'GridAlpha',0.08);
end

title(tl, sprintf(['H2: imposed load at matched ventilation (V_T = %.2f L, RR = %g)   ' ...
    '|  faded+dotted = no published coefficients'], cfg.pattern.V_T_L, cfg.pattern.RR), ...
    'FontWeight','bold','FontSize',12);

viz.save(fig, 'fig_imposed_wob');

% ---------------- H2 verdict ----------------
fprintf('\n  H2 check (matched V_T/RR; grade X devices excluded):\n');
usable = T(~strcmp(T.grade,'X'), :);
for ID = unique(usable.ID_mm)'
    e = usable(strcmp(usable.type,'ETT')   & usable.ID_mm==ID, :);
    t = usable(strcmp(usable.type,'TRACH') & usable.ID_mm==ID, :);
    if ~isempty(e) && ~isempty(t)
        fprintf('    ID %.1f matched: WOB %.3f vs %.3f J/L (trach = %.0f%%) -> %s\n', ...
            ID, e.WOB_J_L, t.WOB_J_L, 100*t.WOB_J_L/e.WOB_J_L, ...
            ternary(t.WOB_J_L < e.WOB_J_L, 'PASS', 'FAIL'));
    end
end

% The clinically realistic pairing is not matched-ID: a trach is usually
% sized up relative to the ETT it replaces.
e75 = T(strcmp(T.device,'ETT_7_5'),:);
t80 = T(strcmp(T.device,'TRACH_8_0'),:);
fprintf('\n    Clinical pairing ETT 7.5 -> TRACH 8.0: WOB %.3f -> %.3f J/L (%.0f%% reduction)\n', ...
    e75.WOB_J_L, t80.WOB_J_L, 100*(1 - t80.WOB_J_L/e75.WOB_J_L));
fprintf('    PTP %.1f -> %.1f cmH2O.s/min\n', e75.PTP_per_min, t80.PTP_per_min);
fprintf('\n    Note: these are IMPOSED loads only. Mode B (exp4) shows what fraction of\n');
fprintf('    TOTAL work this represents -- which is where the difference gets diluted.\n');
end

function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end
