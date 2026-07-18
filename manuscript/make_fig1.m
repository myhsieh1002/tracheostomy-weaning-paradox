function make_fig1()
%MAKE_FIG1 Figure 1 — the four coupled models and the axes they act on.
%
%   Schematic for the manuscript, drawn in the project's viz style and
%   exported to results/figures/fig1_schematic.png so Figure 1 is generated
%   by the same reproducible pipeline as every data figure.
%
%       matlab -batch "addpath('manuscript'); make_fig1"

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
c = viz.style();

fig = figure('Position', [80 80 1180 680], 'Color', 'w');
ax = axes('Position', [0 0 1 1]); hold(ax, 'on');
axis(ax, [0 100 0 100]); axis(ax, 'off');

% ---- palette: one soft tint per model, coloured top rule ----
tint  = @(col) col + (1-col)*0.82;
mB  = c.TRACH;   % Model 1b  — gas exchange (cool)
m1  = c.ETT;     % Model 1   — load (warm)
m2  = c.rescue;  % Model 2   — tipping point (green)
m2b = [0.42 0.33 0.55];  % Model 2b — capacity (violet)

% ---- title / the paradox ----
text(50, 96, 'The tracheostomy weaning–survival paradox', ...
    'HorizontalAlignment','center', 'FontSize',15, 'FontWeight','bold');
text(50, 91.2, ['shortens ventilation and ICU stay  ' char(8212) ...
    '  does not change survival'], ...
    'HorizontalAlignment','center', 'FontSize',10.5, 'Color',c.axisGrey);

% ---- model boxes (axis name is an italic line inside each box) ----
box(4,  56, 22, 22, mB,  'Model 1b', {'CO_2 kinetics','two-compartment'}, 'gas-exchange axis');
box(39, 56, 22, 22, m1,  'Model 1',  {'respiratory mechanics','work of breathing'}, 'load axis');
box(74, 56, 22, 22, m2,  'Model 2',  {'load–capacity','bifurcation'}, 'tipping-point axis');
box(74, 17, 22, 22, m2b, 'Model 2b', {'VIDD capacity','dynamics'}, 'capacity axis');

% ---- couplings ----
arrow(26, 67, 39, 67, c.stable);
lab(32.5, 70.2, {'required V̇_E, V_T'}, c.stable);

arrow(61, 67, 74, 67, c.stable);
lab(67.5, 71.0, {'P_{mus} (load L),', 'f_{device}'}, c.stable);
devtag(67.5, 62.6, {'device lowers load:','small, disease-diluted'}, m1);

arrow(85, 39, 85, 56, c.stable);
lab(80.5, 47.5, {'C_{max}(t)','drives the fold'}, c.stable, 'right');
devtag(89.5, 47.5, {'device may','raise capacity','via sedation','— UNMEASURED'}, m2b, 'left');

% metabolic-rate input into Model 1 (added severity axis)
arrow(50, 84, 50, 78, c.axisGrey);
lab(50, 86, {'metabolic rate (target V̇_A)'}, c.axisGrey);

% ---- bottom thesis strip ----
rectangle('Position',[6 4 88 8.5], 'Curvature',0.5, ...
    'FaceColor',tint(c.grey), 'EdgeColor',c.grey, 'LineWidth',0.75);
text(50, 8.2, ['The device moves the load a little, and least where disease is worst; ' ...
    'it may move capacity a lot, by an amount no one has measured.'], ...
    'HorizontalAlignment','center', 'FontSize',10, 'FontAngle','italic', 'Color',c.stable);

viz.save(fig, 'fig1_schematic');
fprintf('Figure 1 written to results/figures/fig1_schematic.png\n');

% ===================== helpers =====================

    function box(x, y, w, h, col, name, lines, axisCap)
        % soft-tinted rounded box: coloured top rule, title, body, and the
        % axis name as an italic line INSIDE the box (keeps the space below
        % free for coupling labels).
        rectangle('Position',[x y w h], 'Curvature',0.15, ...
            'FaceColor',tint(col), 'EdgeColor',col, 'LineWidth',1.4);
        plot([x+1 x+w-1], [y+h-3.5 y+h-3.5], '-', 'Color',col, 'LineWidth',2.5);
        text(x+w/2, y+h-2.0, name, 'HorizontalAlignment','center', ...
            'VerticalAlignment','middle', 'FontSize',12, 'FontWeight','bold', 'Color',col);
        text(x+w/2, y+h/2+0.5, lines, 'HorizontalAlignment','center', ...
            'VerticalAlignment','middle', 'FontSize',10, 'Color',[0.15 0.15 0.17]);
        text(x+w/2, y+3.0, ['— ' axisCap ' —'], 'HorizontalAlignment','center', ...
            'VerticalAlignment','middle', 'FontSize',9.5, 'FontAngle','italic', 'Color',col);
    end

    function arrow(x1, y1, x2, y2, col)
        plot([x1 x2], [y1 y2], '-', 'Color',col, 'LineWidth',1.6);
        L = hypot(x2-x1, y2-y1); ux = (x2-x1)/L; uy = (y2-y1)/L;
        px = -uy; py = ux; hs = 1.7; hw = 1.0;
        bx = x2 - hs*ux; by = y2 - hs*uy;
        patch([x2 bx+hw*px bx-hw*px], [y2 by+hw*py by-hw*py], col, 'EdgeColor',col);
    end

    function lab(x, y, lines, col, align)
        if nargin < 5, align = 'center'; end
        text(x, y, lines, 'HorizontalAlignment',align, ...
            'VerticalAlignment','middle', 'FontSize',9, 'Color',col);
    end

    function devtag(x, y, txt, col, align)
        if nargin < 5, align = 'center'; end
        text(x, y, txt, 'HorizontalAlignment',align, ...
            'VerticalAlignment','middle', 'FontSize',8.5, 'FontWeight','bold', ...
            'Color',col);
    end
end
