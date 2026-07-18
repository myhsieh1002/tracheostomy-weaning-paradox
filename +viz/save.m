function outPath = save(fig, name)
%SAVE Write a figure to results/figures at publication resolution.
%
%   outPath = viz.save(fig, 'fig_resistance_curves') writes a 300 dpi PNG
%   and a vector PDF alongside it.

arguments
    fig (1,1) matlab.ui.Figure
    name (1,:) char
end

figDir = fullfile(wob.projectRoot, 'results', 'figures');
if ~isfolder(figDir)
    mkdir(figDir);
end

outPath = fullfile(figDir, [name '.png']);
exportgraphics(fig, outPath, 'Resolution', 300);
exportgraphics(fig, fullfile(figDir, [name '.pdf']), 'ContentType', 'vector');
end
