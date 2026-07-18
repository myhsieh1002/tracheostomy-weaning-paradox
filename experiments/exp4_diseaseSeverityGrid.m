function T = exp4_diseaseSeverityGrid(cfg)
%EXP4_DISEASESEVERITYGRID H4 -- the device-fraction dilution. THESIS FIGURE.
%
%   T = exp4_diseaseSeverityGrid() runs Mode B over the (C_rs x R_aw) grid
%   and maps f_device = WOB_device / WOB_total.
%
%   H4: the device's share of the total work of breathing falls
%   monotonically as compliance drops and native resistance rises. In the
%   sickest lungs -- the patients for whom the decision matters most -- the
%   tube is the smallest part of the problem, so removing it moves the least.
%
%   This is the figure the whole Model 1 argument rests on.

arguments
    cfg (1,1) struct = wob.loadConfig()
end

c = viz.style();
Cvals = cfg.disease_grid.C_rs(:)';
Rvals = cfg.disease_grid.R_aw_native(:)';
devices = cellstr(cfg.experiments.exp4.devices);

nC = numel(Cvals); nR = numel(Rvals);
rows = {};
F = struct();
for d = 1:numel(devices)
    F.(devices{d}) = nan(nR, nC);
end

for d = 1:numel(devices)
    for i = 1:nR
        for j = 1:nC
            ov = struct('C_rs', Cvals(j), 'R_aw_native', Rvals(i));
            e = wob.simulateEffort(cfg, devices{d}, ov);
            F.(devices{d})(i,j) = e.f_device;
            rows(end+1,:) = {devices{d}, Cvals(j), Rvals(i), e.f_device, ...
                e.WOB_total_J_min, e.WOB_device_J_min, e.P_mus_mean, e.P_mus_peak, e.V_E}; %#ok<AGROW>
        end
    end
end

T = cell2table(rows, 'VariableNames', ...
    {'device','C_rs','R_aw_native','f_device','WOB_total_J_min','WOB_device_J_min', ...
     'P_mus_mean','P_mus_peak','V_E'});
writetable(T, fullfile(wob.projectRoot,'results','tables','exp4_grid.csv'));

% ---------------- Figure ----------------
fig = figure('Position',[80 80 1180 440]);
tl = tiledlayout(1, numel(devices)+1, 'TileSpacing','compact','Padding','compact');

clim_all = [0, max(cellfun(@(d) max(F.(d)(:)), devices))];

for d = 1:numel(devices)
    nexttile;
    imagesc(Cvals, Rvals, F.(devices{d}));
    set(gca,'YDir','normal'); clim(clim_all);
    colormap(gca, viridisLike());
    xlabel('C_{rs} (L/cmH_2O)'); ylabel('R_{aw,native} (cmH_2O/(L/s))');
    title(viz.deviceLabel(devices{d}));
    xticks(Cvals); yticks(Rvals);

    % Annotate each cell so the figure is readable without a colour bar.
    for i = 1:nR
        for j = 1:nC
            v = F.(devices{d})(i,j);
            text(Cvals(j), Rvals(i), sprintf('%.2f', v), ...
                'HorizontalAlignment','center', 'FontSize', 8, ...
                'Color', ternary(v > 0.5*clim_all(2), [0 0 0], [1 1 1]));
        end
    end
end
cb = colorbar; cb.Label.String = 'f_{device} = WOB_{device} / WOB_{total}';

% --- Dilution curves: f_device vs severity ---
nexttile; hold on;
for d = 1:numel(devices)
    isETT = contains(devices{d},'ETT');
    col = ternary(isETT, c.ETT, c.TRACH);
    for i = 1:nR
        % Lighten towards white for low R_aw, so darker = sicker airway.
        % Interpolating the RGB itself rather than passing an alpha, which
        % plot's Color does not accept.
        mixFrac = 0.65 * (nR - i) / max(nR - 1, 1);
        shade = col + (1 - col) * mixFrac;
        plot(Cvals, F.(devices{d})(i,:), '-o', 'Color', shade, ...
            'MarkerSize',3.5, 'MarkerFaceColor',shade, ...
            'HandleVisibility', ternary(i==nR,'on','off'), ...
            'DisplayName', viz.deviceLabel(devices{d}));
    end
end
xlabel('C_{rs} (L/cmH_2O)  \rightarrow  healthier'); ylabel('f_{device}');
title({'Dilution across severity','(each line one R_{aw}; darker = higher R_{aw})'});
legend('Location','northwest','Box','off');
grid on; set(gca,'GridAlpha',0.08); ylim([0 clim_all(2)*1.05]);

title(tl, 'H4: the device''s share of the work of breathing is diluted by disease', ...
    'FontWeight','bold','FontSize',13);

viz.save(fig, 'fig_device_fraction_heatmap');

% ---------------- H4 monotonicity check ----------------
fprintf('\n  H4 monotonicity check (f_device must fall as disease worsens):\n');
for d = 1:numel(devices)
    M = F.(devices{d});
    % Rows = R_aw ascending; columns = C_rs ascending.
    % f_device must DECREASE as R_aw rises (down columns) and INCREASE as
    % C_rs rises (along rows, i.e. towards a healthier lung).
    monoR = all(diff(M, 1, 1) < 0, 'all');
    monoC = all(diff(M, 1, 2) > 0, 'all');
    fprintf('    %-10s  falls with R_aw: %s | rises with C_rs: %s | range %.3f - %.3f\n', ...
        devices{d}, ternary(monoR,'PASS','FAIL'), ternary(monoC,'PASS','FAIL'), ...
        min(M(:)), max(M(:)));
end
end

function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end

function m = viridisLike()
% A perceptually ordered map without requiring a toolbox colormap.
anchors = [0.267 0.005 0.329; 0.283 0.141 0.458; 0.254 0.265 0.530;
           0.207 0.372 0.553; 0.164 0.471 0.558; 0.128 0.567 0.551;
           0.135 0.659 0.518; 0.267 0.749 0.441; 0.478 0.821 0.318;
           0.741 0.873 0.150; 0.993 0.906 0.144];
x = linspace(0,1,size(anchors,1));
xi = linspace(0,1,256);
m = [interp1(x,anchors(:,1),xi)', interp1(x,anchors(:,2),xi)', interp1(x,anchors(:,3),xi)'];
end
