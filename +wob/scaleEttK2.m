function cfg = scaleEttK2(cfg, k)
%SCALEETTK2 Represent partially obstructed ("in-use") endotracheal tubes.
%
%   cfg = wob.scaleEttK2(cfg, k) multiplies the quadratic Rohrer coefficient
%   K2 of every ETT in the ACTIVE coefficient set by k. Tracheostomy tubes
%   and the native upper airway are left unchanged.
%
%   Bench coefficients describe new tubes. Flevari 2011 (PMID 21675060)
%   measured in-vivo ETT K2 1.3-4.3x the bench value after 4-11 days of
%   intubation, and Diehl 1999 (PMID 9927347) found used ETTs carried 28%
%   more resistive work than new ones. k = 1 reproduces the primary model.
%   K1 is left alone: it is ~0 for the active set and Flevari's in-vivo
%   excess is in the quadratic term.
%
%   See also wob.loadConfig, m2_exp6_inuseETT

arguments
    cfg (1,1) struct
    k   (1,1) double {mustBePositive}
end

act = cfg.options_coefficient_set.active;
names = fieldnames(cfg.devices);
for i = 1:numel(names)
    d = cfg.devices.(names{i});
    if strcmp(d.type, 'ETT')
        cfg.devices.(names{i}).coefficients.(act).K2 = d.coefficients.(act).K2 * k;
    end
end
end
