function c = style()
%STYLE Shared figure styling and the project colour map.
%
%   c = viz.style() returns a struct of colours and applies publication
%   defaults to the current graphics root.
%
%   Device colours are held fixed across every figure so that ETT and
%   tracheostomy read the same way in each panel of the final triptych.

set(groot, 'DefaultAxesFontSize', 11);
set(groot, 'DefaultAxesLineWidth', 0.75);
set(groot, 'DefaultLineLineWidth', 1.6);
set(groot, 'DefaultAxesBox', 'off');
set(groot, 'DefaultAxesTickDir', 'out');
set(groot, 'DefaultFigureColor', 'w');

% Device identity: warm = ETT (the incumbent), cool = tracheostomy.
c.ETT        = [0.84 0.37 0.20];
c.ETT_light  = [0.94 0.68 0.55];
c.TRACH      = [0.16 0.44 0.68];
c.TRACH_light= [0.55 0.74 0.88];

% Work decomposition
c.elastic    = [0.45 0.45 0.48];
c.native     = [0.62 0.72 0.45];
c.device     = [0.84 0.37 0.20];

% Bifurcation structure
c.stable     = [0.13 0.13 0.15];
c.unstable   = [0.70 0.70 0.72];
c.fold       = [0.72 0.15 0.25];
c.rescue     = [0.22 0.55 0.35];

c.grey       = [0.55 0.55 0.58];
c.axisGrey   = [0.35 0.35 0.38];
end
