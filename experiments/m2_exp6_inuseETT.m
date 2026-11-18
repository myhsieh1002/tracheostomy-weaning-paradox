function S = m2_exp6_inuseETT(cfg1, cfg2)
%M2_EXP6_INUSEETT Sensitivity: what if the endotracheal tube is not new?
%
%   S = m2_exp6_inuseETT() repeats the two headline computations of Paper A
%   with every ETT's quadratic coefficient scaled by k (wob.scaleEttK2):
%
%     - the device fraction f_device over the exp4 disease grid (ETT 7.5),
%       and whether its monotone dilution with disease survives;
%     - the rescue-window width dL_device / l_high over the m2_exp3 grid
%       (R_aw = 10, three metabolic levels), and its overlap with the
%       observed weaning capacity range (Vassilakopoulos MIP 42.3-53.8);
%     - the fall in TOTAL work of breathing from ETT 7.5 to tracheostomy 8.0,
%       for the comparison against Diehl 1999's within-patient -55%.
%
%   WHY
%   ---
%   The primary model uses bench coefficients of NEW tubes. Flevari 2011
%   measured in-vivo ETT K2 1.3-4.3x the bench value after 4-11 days of
%   intubation -- the window in which tracheostomy is decided -- which is
%   beyond the Sobol K2_scale range of exp5 (upper bound 1.5). This makes
%   the direction of that bias explicit: the primary device fraction and
%   rescue-band width are lower bounds for a bare tracheostomy tube. (An
%   inner cannula biases the other way -- Carter 2013, Pryor 2016 -- and is
%   not modelled.)
%
%   See also wob.scaleEttK2, exp4_diseaseSeverityGrid, m2_exp3_rescueWindow

arguments
    cfg1 (1,1) struct = wob.loadConfig()
    cfg2 (1,1) struct = lc.loadConfig()
end

ks = reshape(cfg1.inuse_ett.K2_scale, 1, []);
p  = lc.calibrate(cfg2);
nf = lc.normalizedFolds(p);
devETT = cfg2.load_from_model1.devices{1};
devTR  = cfg2.load_from_model1.devices{2};
MIP    = [cfg2.validation.vassilakopoulos1998.failure.MIP_cmH2O, ...
          cfg2.validation.vassilakopoulos1998.success.MIP_cmH2O];

Cv = cfg1.disease_grid.C_rs(:)';
Rv = cfg1.disease_grid.R_aw_native(:)';
VA = cfg1.disease_grid.target_VA_L_min(:)';
Cvals   = linspace(0.02, 0.06, 110);      % as m2_exp3
R_fixed = 10;                             % as m2_exp3

summ = zeros(numel(ks), 12);
byVA = [];
for ik = 1:numel(ks)
    k = ks(ik);
    c1 = wob.scaleEttK2(cfg1, k);

    % ---- Model 1: f_device and total-WOB fall over the exp4 grid ----
    F = nan(numel(Rv), numel(Cv)); D = F;
    for i = 1:numel(Rv)
        for j = 1:numel(Cv)
            o  = struct('C_rs', Cv(j), 'R_aw_native', Rv(i));
            eE = wob.simulateEffort(c1, 'ETT_7_5', o);
            eT = wob.simulateEffort(c1, 'TRACH_8_0', o);
            F(i,j) = eE.f_device;
            D(i,j) = 1 - eT.WOB_total_J_min / eE.WOB_total_J_min;
        end
    end
    monoR = all(diff(F,1,1) < 0, 'all');   % falls as airway resistance rises
    monoC = all(diff(F,1,2) > 0, 'all');   % rises as compliance rises (healthier)

    % ---- Model 2: rescue window, as m2_exp3 ----
    W = [];
    for v = 1:numel(VA)
        cv = c1; cv.patient.target_VA_L_min = VA(v);
        LE = zeros(size(Cvals)); LT = LE;
        for j = 1:numel(Cvals)
            o = struct('C_rs', Cvals(j), 'R_aw_native', R_fixed);
            LE(j) = lc.coupling(cv, cfg2, devETT, o).L_total;
            LT(j) = lc.coupling(cv, cfg2, devTR,  o).L_total;
        end
        w  = (LE - LT) / nf.l_high;
        ov = any(LE/nf.l_high > MIP(1) & LT/nf.l_high < MIP(2));
        W  = [W, w]; %#ok<AGROW>
        byVA = [byVA; k, VA(v), min(w), max(w), mean(w), ov]; %#ok<AGROW>
    end

    summ(ik,:) = [k, min(F(:)), max(F(:)), monoR && monoC, ...
                  mean(W), min(W), max(W), mean(W)/50, ...
                  D(1,end), D(end,1), min(D(:)), max(D(:))];
end

S.summary = array2table(summ, 'VariableNames', {'K2_scale','f_device_min','f_device_max', ...
    'dilution_monotone','window_mean_cmH2O','window_min_cmH2O','window_max_cmH2O', ...
    'window_frac_of_20_70','WOB_fall_healthiest','WOB_fall_sickest','WOB_fall_min','WOB_fall_max'});
S.byVA = array2table(byVA, 'VariableNames', {'K2_scale','target_VA','window_min_cmH2O', ...
    'window_max_cmH2O','window_mean_cmH2O','overlaps_observed_MIP'});

outDir = fullfile(wob.projectRoot, 'results', 'tables');
writetable(S.summary, fullfile(outDir, 'm2_exp6_inuseETT_summary.csv'));
writetable(S.byVA,    fullfile(outDir, 'm2_exp6_inuseETT_window_by_VA.csv'));

fprintf('\n  In-use ETT scenario (ETT K2 x k; tracheostomy unchanged):\n');
disp(S.summary);
end
