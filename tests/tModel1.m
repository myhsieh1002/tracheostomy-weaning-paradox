classdef tModel1 < matlab.unittest.TestCase
%TMODEL1 Model 1: devices, waveforms, mechanics, metrics.
%
%   Covers the sanity checks in section 7 of the Model 1 build plan, each
%   against a closed form or a published number rather than against the
%   code's own previous output.

    properties
        cfg
    end

    methods (TestClassSetup)
        function setup(t)
            addpath(fileparts(fileparts(mfilename('fullpath'))));
            t.cfg = wob.loadConfig();
        end
    end

    methods (Test)

        % ---------------- units ----------------

        function unitConversionIsExact(t)
            % Spec section 7 asks this to be hard-checked, because getting it
            % wrong would rescale every reported joule silently.
            c = wob.constants();
            t.verifyEqual(c.CMH2O_L_TO_JOULE, 0.0981, 'AbsTol', 1e-12);
            t.verifyEqual(c.MMHG_TO_CMH2O, 1.35951, 'AbsTol', 1e-5);
        end

        % ---------------- devices ----------------

        function rohrerIsMonotoneInFlow(t)
            flow = linspace(0, 2, 200);
            names = setdiff(fieldnames(t.cfg.devices), {'NATIVE_UPPER_AIRWAY'});
            for k = 1:numel(names)
                dev = wob.getDevice(t.cfg, names{k});
                dP = wob.rohrerDrop(dev, flow);
                t.verifyTrue(all(diff(dP) > -1e-12), ...
                    sprintf('%s: pressure drop must be non-decreasing in flow', names{k}));
            end
        end

        function rohrerIsOddSymmetric(t)
            % Expiratory flow must give an equal and opposite drop, or the
            % expiratory limb would be modelled with the wrong sign.
            dev = wob.getDevice(t.cfg, 'ETT_7_5');
            t.verifyEqual(wob.rohrerDrop(dev, -0.8), -wob.rohrerDrop(dev, 0.8), 'RelTol', 1e-12);
        end

        function trachHasLowerResistanceAtMatchedID(t)
            % H1. Only ID 8.0 has published coefficients for both classes.
            e = wob.getDevice(t.cfg, 'ETT_8_0');
            r = wob.getDevice(t.cfg, 'TRACH_8_0');
            t.verifyLessThan(wob.rohrerDrop(r, 1.0), wob.rohrerDrop(e, 1.0));
        end

        function guttmannWithinRigRatioReproduced(t)
            % Guttmann 1993 reports the 8.0 trach at ~70% of the 8.0 ETT at
            % 1 L/s. Reproducing his own published check validates the refit
            % and the transcription end to end.
            t.assumeEqual(t.cfg.options_coefficient_set.active, 'guttmann_within_rig');
            e = wob.getDevice(t.cfg, 'ETT_8_0');
            r = wob.getDevice(t.cfg, 'TRACH_8_0');
            ratio = wob.rohrerDrop(r, 1.0) / wob.rohrerDrop(e, 1.0);
            t.verifyEqual(ratio, 0.70, 'AbsTol', 0.03);
        end

        function calibrateRohrerRecoversKnownCoefficients(t)
            % Round-trip: synthesise from known K1/K2, fit, recover them.
            K1true = 1.7; K2true = 8.3;
            flow = linspace(0.1, 1.5, 25)';
            dP = K1true*flow + K2true*flow.*abs(flow);
            [K1, K2, info] = wob.calibrateRohrer(flow, dP);
            t.verifyEqual(K1, K1true, 'RelTol', 1e-10);
            t.verifyEqual(K2, K2true, 'RelTol', 1e-10);
            t.verifyGreaterThan(info.R2, 1 - 1e-12);
        end

        function calibrateRohrerRejectsMismatchedInput(t)
            t.verifyError(@() wob.calibrateRohrer([1;2;3], [1;2]), ...
                'wob:calibrateRohrer:sizeMismatch');
        end

        function physicsScalingExponentIsInTurbulentRange(t)
            % Spec section 7: resistance should scale as roughly 1/d^4..5.
            IDs = linspace(6, 9, 40);
            dP = arrayfun(@(d) abs(wob.physicsScaling(d, 30, 1.0).dP_total), IDs);
            expo = mean(diff(log(dP)) ./ diff(log(IDs)));
            t.verifyGreaterThan(expo, -5.0);
            t.verifyLessThan(expo, -4.0);
        end

        % ---------------- waveform ----------------

        function inspiratoryFlowIntegratesToTidalVolume(t)
            for wf = ["sinusoidal","constant","decelerating","accelerating"]
                pat = t.cfg.pattern; pat.waveform = wf;
                b = wob.breathingPattern(pat);
                V = trapz(b.t(b.isInsp), b.flow(b.isInsp));
                t.verifyEqual(V, pat.V_T_L, 'RelTol', 1e-12, ...
                    sprintf('waveform %s must deliver exactly V_T', wf));
            end
        end

        function halfSineIsAnAliasOfSinusoidal(t)
            % Documented as an alias rather than invented as a distinct
            % shape; assert that it stays one.
            pa = t.cfg.pattern; pa.waveform = 'sinusoidal';
            pb = t.cfg.pattern; pb.waveform = 'half_sine';
            t.verifyEqual(wob.breathingPattern(pb).flow, wob.breathingPattern(pa).flow, 'RelTol', 1e-12);
        end

        function dutyCycleFollowsIERatio(t)
            pat = t.cfg.pattern; pat.IE_ratio = [1 2];
            b = wob.breathingPattern(pat);
            t.verifyEqual(b.dutyCycle, 1/3, 'RelTol', 1e-9);
            pat.IE_ratio = [1 3];
            t.verifyEqual(wob.breathingPattern(pat).dutyCycle, 0.25, 'RelTol', 1e-9);
        end

        function unknownWaveformErrors(t)
            pat = t.cfg.pattern; pat.waveform = 'triangular';
            t.verifyError(@() wob.breathingPattern(pat), 'wob:breathingPattern:unknownWaveform');
        end

        % ---------------- dead space ----------------

        function alveolarVentilationIdentityHolds(t)
            % V_A = (V_T - V_D) * RR, numerically.
            dev = wob.getDevice(t.cfg, 'TRACH_8_0');
            v = wob.requiredVentilation(t.cfg, dev);
            t.verifyEqual((v.V_T - v.V_D) * v.RR, t.cfg.patient.target_VA_L_min, 'RelTol', 1e-12);
        end

        function bypassTermCancelsBetweenETTandTrach(t)
            % The central dead-space finding: both devices bypass the same
            % region, so dV_D cannot influence the ETT-vs-trach contrast.
            % Changing it must move both dead spaces equally.
            cfgA = t.cfg; cfgA.patient.Vd_upper_bypassed_mL = 40;
            cfgB = t.cfg; cfgB.patient.Vd_upper_bypassed_mL = 90;
            dA = wob.deadSpace(cfgA, wob.getDevice(cfgA,'ETT_7_5')).total ...
               - wob.deadSpace(cfgA, wob.getDevice(cfgA,'TRACH_8_0')).total;
            dB = wob.deadSpace(cfgB, wob.getDevice(cfgB,'ETT_7_5')).total ...
               - wob.deadSpace(cfgB, wob.getDevice(cfgB,'TRACH_8_0')).total;
            t.verifyEqual(dA, dB, 'AbsTol', 1e-12, ...
                'dV_D must cancel from the ETT-vs-trach dead-space difference');
        end

        function deviceDeadSpaceDifferenceIsApparatusOnly(t)
            e = wob.getDevice(t.cfg,'ETT_7_5');
            r = wob.getDevice(t.cfg,'TRACH_8_0');
            diffTotal = wob.deadSpace(t.cfg, e).total - wob.deadSpace(t.cfg, r).total;
            diffApp = (e.Vd_apparatus_mL - r.Vd_apparatus_mL) * 1e-3;
            t.verifyEqual(diffTotal, diffApp, 'AbsTol', 1e-12);
        end

        function bypassExceedingAnatomicErrors(t)
            cfg2 = t.cfg; cfg2.patient.Vd_upper_bypassed_mL = 500;
            t.verifyError(@() wob.deadSpace(cfg2, wob.getDevice(cfg2,'ETT_7_5')), ...
                'wob:deadSpace:bypassExceedsAnatomic');
        end

        % ---------------- mechanics & metrics ----------------

        function elasticWorkMatchesClosedForm(t)
            % WOB_elastic must equal 0.5*V_T^2/C exactly.
            dev = wob.getDevice(t.cfg, 'ETT_7_5');
            r = wob.simulateModeB(t.cfg, dev);
            c = wob.constants();
            expected = 0.5 * r.V_T^2 / t.cfg.patient.C_rs * c.CMH2O_L_TO_JOULE;
            t.verifyEqual(r.WOB_elastic_J, expected, 'RelTol', 1e-9);
        end

        function workDecompositionSumsToTotal(t)
            dev = wob.getDevice(t.cfg, 'ETT_7_5');
            r = wob.simulateModeB(t.cfg, dev);
            parts = r.WOB_elastic_J + r.WOB_native_J + r.WOB_device_J + ...
                    r.WOB_inertial_J + r.WOB_offset_J;
            t.verifyEqual(parts, r.WOB_total_J, 'RelTol', 1e-12);
        end

        function resistiveWorkMatchesClosedFormForConstantFlow(t)
            % With square flow and a purely linear device, the resistive work
            % is analytic: W = R*Vdot*V_T with Vdot = V_T/Ti.
            cfg2 = t.cfg;
            cfg2.pattern.waveform = 'constant';
            cfg2.patient.R_aw_native = 5;
            activeSet = cfg2.options_coefficient_set.active;
            cfg2.devices.ETT_7_5.coefficients.(activeSet).K1 = 0;
            cfg2.devices.ETT_7_5.coefficients.(activeSet).K2 = 0;
            dev = wob.getDevice(cfg2, 'ETT_7_5');
            r = wob.simulateModeB(cfg2, dev);
            c = wob.constants();
            Vdot = r.V_T / r.pat.Ti;
            expected = 5 * Vdot * r.V_T * c.CMH2O_L_TO_JOULE;
            t.verifyEqual(r.WOB_native_J, expected, 'RelTol', 1e-9);
        end

        function deviceWorkVanishesForZeroResistance(t)
            % Spec section 7 limit: as the tube stops resisting,
            % WOB_device -> 0 and f_device -> 0.
            cfg2 = t.cfg;
            activeSet = cfg2.options_coefficient_set.active;
            cfg2.devices.ETT_7_5.coefficients.(activeSet).K1 = 0;
            cfg2.devices.ETT_7_5.coefficients.(activeSet).K2 = 0;
            e = wob.simulateEffort(cfg2, 'ETT_7_5');
            t.verifyEqual(e.WOB_device, 0, 'AbsTol', 1e-12);
            t.verifyEqual(e.f_device, 0, 'AbsTol', 1e-12);
        end

        function fDeviceRisesWithDeviceResistance(t)
            cfg2 = t.cfg;
            activeSet = cfg2.options_coefficient_set.active;
            base = wob.simulateEffort(cfg2, 'ETT_7_5').f_device;
            cfg2.devices.ETT_7_5.coefficients.(activeSet).K2 = ...
                cfg2.devices.ETT_7_5.coefficients.(activeSet).K2 * 3;
            t.verifyGreaterThan(wob.simulateEffort(cfg2,'ETT_7_5').f_device, base);
        end

        function ptpIsPositiveAndScalesWithRate(t)
            dev = wob.getDevice(t.cfg, 'ETT_7_5');
            r = wob.simulateModeB(t.cfg, dev);
            t.verifyGreaterThan(r.PTP_per_breath, 0);
            t.verifyEqual(r.PTP_per_min, r.PTP_per_breath * r.RR, 'RelTol', 1e-12);
        end

        function pmusSummariesAreOrdered(t)
            % peak >= mean >= duty-weighted, by construction.
            dev = wob.getDevice(t.cfg, 'ETT_7_5');
            r = wob.simulateModeB(t.cfg, dev);
            t.verifyGreaterThanOrEqual(r.P_mus_peak, r.P_mus_mean);
            t.verifyGreaterThanOrEqual(r.P_mus_mean, r.P_mus_dutyWeighted);
            t.verifyEqual(r.P_mus_dutyWeighted, r.P_mus_mean * r.pat.dutyCycle, 'RelTol', 1e-9);
        end

        % ---------------- H4 ----------------

        function fDeviceFallsMonotonicallyWithDiseaseSeverity(t)
            % H4, asserted over the whole configured grid rather than at a
            % couple of points.
            Cvals = t.cfg.disease_grid.C_rs(:)';
            Rvals = t.cfg.disease_grid.R_aw_native(:)';
            F = zeros(numel(Rvals), numel(Cvals));
            for i = 1:numel(Rvals)
                for j = 1:numel(Cvals)
                    F(i,j) = wob.simulateEffort(t.cfg, 'ETT_7_5', ...
                        struct('C_rs',Cvals(j),'R_aw_native',Rvals(i))).f_device;
                end
            end
            t.verifyTrue(all(diff(F,1,1) < 0, 'all'), 'f_device must fall as R_aw rises');
            t.verifyTrue(all(diff(F,1,2) > 0, 'all'), 'f_device must rise as C_rs rises');
        end

        % ---------------- config integrity ----------------

        function unknownOverrideIsRejected(t)
            t.verifyError(@() wob.simulateEffort(t.cfg,'ETT_7_5',struct('C_rz',0.03)), ...
                'wob:simulateEffort:unknownOverride');
        end

        function unknownDeviceIsRejected(t)
            t.verifyError(@() wob.getDevice(t.cfg,'ETT_9_9'), 'wob:getDevice:unknownDevice');
        end

        function unknownCoefficientSetIsRejected(t)
            t.verifyError(@() wob.getDevice(t.cfg,'ETT_7_5','nonesuch'), ...
                'wob:getDevice:unknownCoefficientSet');
        end

        function everyDeviceCarriesProvenance(t)
            names = fieldnames(t.cfg.devices);
            for k = 1:numel(names)
                dev = wob.getDevice(t.cfg, names{k});
                t.verifyNotEmpty(dev.source, sprintf('%s has no source', names{k}));
                t.verifyTrue(ismember(dev.grade, {'A','B','X','A/B'}), ...
                    sprintf('%s has an invalid evidence grade "%s"', names{k}, dev.grade));
            end
        end

        function conclusionIsRobustToCoefficientSet(t)
            % Swapping the whole bench source must not change f_device much;
            % if it did, the result would be an artefact of the source rather
            % than of the physiology.
            a = t.cfg; a.options_coefficient_set.active = 'guttmann_within_rig';
            b = t.cfg; b.options_coefficient_set.active = 'flevari_ett';
            fa = wob.simulateEffort(a,'ETT_7_5').f_device;
            fb = wob.simulateEffort(b,'ETT_7_5').f_device;
            t.verifyEqual(fa, fb, 'AbsTol', 0.05);
        end
    end
end
