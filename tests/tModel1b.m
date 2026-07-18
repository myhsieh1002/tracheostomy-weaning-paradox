classdef tModel1b < matlab.unittest.TestCase
%TMODEL1B Model 1b: CO2 kinetics, ventilatory demand, coupling.
%
%   Covers section 6 of the Model 1b build plan, plus the two structural
%   claims this model rests on: that dead space is single-sourced with
%   Model 1, and that the alveolar dead space is parametrised so the device
%   does not vanish from the answer.

    properties
        cfgCO2
        cfgM1
        cfgM2
    end

    methods (TestClassSetup)
        function setup(t)
            addpath(fileparts(fileparts(mfilename('fullpath'))));
            t.cfgCO2 = co2.loadConfig();
            t.cfgM1  = wob.loadConfig();
            t.cfgM2  = lc.loadConfig();
        end
    end

    methods (Test)

        % ---------------- the alveolar gas equation ----------------

        function alveolarConstantReproducesTargetVA(t)
            % The constant that ties this model to Model 1: VCO2 = 200 and
            % PaCO2 = 40 must give V_A ~ 4.3 L/min, which is where Model 1's
            % target_VA came from.
            K = t.cfgCO2.constants.alveolar_K;
            V_A = K * 200 / 40;
            t.verifyEqual(V_A, 4.315, 'AbsTol', 0.02);
        end

        function steadyStateInvertsRequiredVentilation(t)
            % The two directions must be exact inverses of each other.
            v = co2.requiredVentilation(t.cfgCO2, t.cfgM1, 'ETT_7_5');
            s = co2.steadyState(t.cfgCO2, t.cfgM1, 'ETT_7_5', V_T=v.V_T, RR=v.RR);
            t.verifyEqual(s.PaCO2, t.cfgCO2.targets.PaCO2_target_mmHg, 'RelTol', 1e-10);
        end

        function paCO2ScalesWithVCO2AtFixedVentilation(t)
            % PaCO2 proportional to VCO2 -- spec section 6 sanity check.
            a = co2.steadyState(t.cfgCO2, t.cfgM1, 'ETT_7_5', VCO2=200);
            b = co2.steadyState(t.cfgCO2, t.cfgM1, 'ETT_7_5', VCO2=400);
            t.verifyEqual(b.PaCO2, 2*a.PaCO2, 'RelTol', 1e-12);
        end

        function paCO2ScalesInverselyWithAlveolarVentilation(t)
            a = co2.steadyState(t.cfgCO2, t.cfgM1, 'ETT_7_5', RR=18);
            b = co2.steadyState(t.cfgCO2, t.cfgM1, 'ETT_7_5', RR=36);
            t.verifyEqual(a.PaCO2 * a.V_A, b.PaCO2 * b.V_A, 'RelTol', 1e-12);
        end

        function requiredVentilationScalesWithVCO2(t)
            a = co2.requiredVentilation(t.cfgCO2, t.cfgM1, 'ETT_7_5', VCO2=200);
            b = co2.requiredVentilation(t.cfgCO2, t.cfgM1, 'ETT_7_5', VCO2=400);
            t.verifyEqual(b.V_A, 2*a.V_A, 'RelTol', 1e-12);
        end

        % ---------------- dead space: single source, and the device survives ----------------

        function seriesDeadSpaceComesFromModel1(t)
            % Spec section 8: dV_D must be single-sourced. Assert this model
            % returns Model 1's number rather than its own.
            dev = wob.getDevice(t.cfgM1, 'ETT_7_5');
            m1  = wob.deadSpace(t.cfgM1, dev);
            ds  = co2.deadSpaceTotal(t.cfgCO2, t.cfgM1, 'ETT_7_5', 0.45);
            t.verifyEqual(ds.series, m1.total, 'AbsTol', 1e-15);
        end

        function changingModel1BypassPropagatesHere(t)
            % The single-source claim has to survive a change at the source.
            cfgA = t.cfgM1; cfgA.patient.Vd_upper_bypassed_mL = 40;
            cfgB = t.cfgM1; cfgB.patient.Vd_upper_bypassed_mL = 90;
            dsA = co2.deadSpaceTotal(t.cfgCO2, cfgA, 'ETT_7_5', 0.45);
            dsB = co2.deadSpaceTotal(t.cfgCO2, cfgB, 'ETT_7_5', 0.45);
            t.verifyNotEqual(dsA.series, dsB.series);
            t.verifyEqual(dsA.series - dsB.series, 50e-3, 'AbsTol', 1e-12);
        end

        function bypassStillCancelsBetweenETTandTrach(t)
            % Model 1's central dead-space finding must be inherited, not
            % undone by the alveolar term: dV_D moves both arms equally, so
            % the ETT-vs-trach ventilatory difference cannot depend on it.
            phi = co2.alveolarDeadSpaceFraction(t.cfgCO2, t.cfgM1);
            bypasses = [40 90];
            d = zeros(1, numel(bypasses));
            for k = 1:numel(bypasses)
                c1 = t.cfgM1;
                c1.patient.Vd_upper_bypassed_mL = bypasses(k);
                vE = co2.requiredVentilation(t.cfgCO2, c1, 'ETT_7_5',   Phi=phi);
                vT = co2.requiredVentilation(t.cfgCO2, c1, 'TRACH_8_0', Phi=phi);
                d(k) = vE.V_E - vT.V_E;
            end
            t.verifyEqual(d(1), d(2), 'RelTol', 1e-12, ...
                'dV_D must cancel from the ETT-vs-trach ventilatory difference');
        end

        function deviceStillAffectsVentilatoryDemand(t)
            % THE regression test for this model's central modelling fix.
            % If the alveolar dead space were made a fraction of V_T rather
            % than of alveolar gas, V_E = V_A_req/(1-VD_VT) would contain no
            % device term at all and this would return zero.
            vE = co2.requiredVentilation(t.cfgCO2, t.cfgM1, 'ETT_7_5');
            vT = co2.requiredVentilation(t.cfgCO2, t.cfgM1, 'TRACH_8_0');
            t.verifyGreaterThan(vE.V_E - vT.V_E, 1e-6, ...
                'the device must not vanish from ventilatory demand');
        end

        function deviceEffectEqualsApparatusDifferenceTimesRate(t)
            % dV_E = (V_D_series,ETT - V_D_series,trach) * RR, exactly.
            phi = co2.alveolarDeadSpaceFraction(t.cfgCO2, t.cfgM1);
            vE = co2.requiredVentilation(t.cfgCO2, t.cfgM1, 'ETT_7_5',   Phi=phi);
            vT = co2.requiredVentilation(t.cfgCO2, t.cfgM1, 'TRACH_8_0', Phi=phi);
            dSeries = vE.ds.series - vT.ds.series;
            t.verifyEqual(vE.V_E - vT.V_E, dSeries * vE.RR, 'RelTol', 1e-9);
        end

        function alveolarFractionReproducesConfiguredVDVT(t)
            % phi is calibrated so the configured V_D/V_T comes back at the
            % reference condition.
            phi = co2.alveolarDeadSpaceFraction(t.cfgCO2, t.cfgM1);
            ds  = co2.deadSpaceTotal(t.cfgCO2, t.cfgM1, t.cfgCO2.arms.ett, ...
                                     t.cfgCO2.pattern.V_T_L, phi);
            t.verifyEqual(ds.VD_VT_eff, t.cfgCO2.deadspace.VD_VT, 'RelTol', 1e-9);
        end

        function impossibleVDVTIsRejected(t)
            % A V_D/V_T below what the series dead space alone implies cannot
            % be achieved by any V/Q distribution; reject rather than clamp.
            c = t.cfgCO2; c.deadspace.VD_VT = 0.05;
            t.verifyError(@() co2.alveolarDeadSpaceFraction(c, t.cfgM1), ...
                'co2:alveolarDeadSpaceFraction:impossibleVDVT');
        end

        function unityVDVTIsRejected(t)
            c = t.cfgCO2; c.deadspace.VD_VT = 1.0;
            t.verifyError(@() co2.alveolarDeadSpaceFraction(c, t.cfgM1), ...
                'co2:alveolarDeadSpaceFraction:noAlveolarVentilation');
        end

        function zeroDeadSpaceLimitGivesVEequalsVA(t)
            % Spec section 6: as V_D -> 0, V_E -> V_A.
            c1 = t.cfgM1;
            c1.patient.Vd_anatomic_mL = 0;
            c1.patient.Vd_upper_bypassed_mL = 0;
            c1.devices.ETT_7_5.Vd_apparatus_mL = 0;
            c = t.cfgCO2;
            c.deadspace.VD_VT = 1e-9;   % phi -> 0
            v = co2.requiredVentilation(c, c1, 'ETT_7_5');
            t.verifyEqual(v.V_E, v.V_A, 'RelTol', 1e-6);
        end

        % ---------------- H3 dilution ----------------

        function benefitShareFallsWithMetabolicRate(t)
            phi = co2.alveolarDeadSpaceFraction(t.cfgCO2, t.cfgM1);
            frac = zeros(1,2); VCO2s = [200 400];
            for k = 1:2
                vE = co2.requiredVentilation(t.cfgCO2, t.cfgM1, 'ETT_7_5',   VCO2=VCO2s(k), Phi=phi);
                vT = co2.requiredVentilation(t.cfgCO2, t.cfgM1, 'TRACH_8_0', VCO2=VCO2s(k), Phi=phi);
                frac(k) = (vE.V_E - vT.V_E) / vE.V_E;
            end
            t.verifyLessThan(frac(2), frac(1), 'the dead-space share must dilute as VCO2 rises');
        end

        function benefitShareFallsWithVQMismatch(t)
            frac = zeros(1,2); VDVTs = [0.3 0.6];
            for k = 1:2
                c = t.cfgCO2; c.deadspace.VD_VT = VDVTs(k);
                phi = co2.alveolarDeadSpaceFraction(c, t.cfgM1);
                vE = co2.requiredVentilation(c, t.cfgM1, 'ETT_7_5',   Phi=phi);
                vT = co2.requiredVentilation(c, t.cfgM1, 'TRACH_8_0', Phi=phi);
                frac(k) = (vE.V_E - vT.V_E) / vE.V_E;
            end
            t.verifyLessThan(frac(2), frac(1), 'the dead-space share must dilute as V/Q worsens');
        end

        function bypassContrastIsMuchLargerThanDeviceContrast(t)
            % The 72 mL literature is about the first, clinicians choose the
            % second. They must not be conflated.
            phi = co2.alveolarDeadSpaceFraction(t.cfgCO2, t.cfgM1);
            vN = co2.requiredVentilation(t.cfgCO2, t.cfgM1, 'NATIVE_UPPER_AIRWAY', Phi=phi);
            vE = co2.requiredVentilation(t.cfgCO2, t.cfgM1, 'ETT_7_5',             Phi=phi);
            vT = co2.requiredVentilation(t.cfgCO2, t.cfgM1, 'TRACH_8_0',           Phi=phi);
            t.verifyGreaterThan(vN.V_E - vT.V_E, 3 * (vE.V_E - vT.V_E));
        end

        % ---------------- two-compartment transient ----------------

        function transientSettlesToTheAlgebraicSteadyState(t)
            % The ODE's steady state must reproduce the alveolar gas
            % equation, or the unit chain is wrong somewhere.
            V_A = 4.315;
            out = co2.twoCompartment(t.cfgCO2, [0 600], @(tt) V_A, VCO2=200);
            t.verifyEqual(out.PA(end), 40, 'RelTol', 5e-3);
        end

        function venousArterialDifferenceIsPhysiologic(t)
            % Pv - PA = VCO2/(Q*slope*1000); ~6 mmHg is the physiologic value
            % and is the check that the dissociation slope is the right order.
            out = co2.twoCompartment(t.cfgCO2, [0 600], @(tt) 4.315, VCO2=200);
            t.verifyEqual(out.Pv(end) - out.PA(end), 6.67, 'AbsTol', 0.5);
        end

        function bodyTimeConstantIsMinutesNotSeconds(t)
            % The spec's bare S_body = 15 read as an absolute capacitance
            % gives a 0.5 min body time constant, contradicting its own
            % section 2.1. Reinterpreted per-kg it should be tens of minutes.
            out = co2.twoCompartment(t.cfgCO2, [0 10], @(tt) 4.315);
            t.verifyGreaterThan(out.tau_body, 10, 'body CO2 stores must be slow (tens of minutes)');
            t.verifyLessThan(out.tau_body, 120);
        end

        function lungTimeConstantIsSeconds(t)
            out = co2.twoCompartment(t.cfgCO2, [0 10], @(tt) 4.315);
            t.verifyLessThan(out.tau_lung, 2, 'lung washout must be fast (under ~1 min)');
            t.verifyGreaterThan(out.tau_lung, 0.1);
        end

        function timescalesAreSeparated(t)
            % The whole justification for two compartments.
            out = co2.twoCompartment(t.cfgCO2, [0 10], @(tt) 4.315);
            t.verifyGreaterThan(out.tau_body / out.tau_lung, 10);
        end

        function raisingVentilationLowersPaCO2Transiently(t)
            out = co2.twoCompartment(t.cfgCO2, [0 240], @(tt) 4.315 * (1 + 0.5*(tt > 30)));
            iBefore = find(out.t < 30, 1, 'last');
            t.verifyLessThan(out.PA(end), out.PA(iBefore));
        end

        % ---------------- coupling ----------------

        function couplingFeedsModel1TheVentilationItComputed(t)
            cp = co2.coupling(t.cfgCO2, t.cfgM1, t.cfgM2, 'ETT_7_5');
            t.verifyEqual(cp.effort.V_T, cp.vent.V_T, 'RelTol', 1e-12);
            t.verifyEqual(cp.effort.RR,  cp.vent.RR,  'RelTol', 1e-12);
        end

        function model1DoesNotResolveTheVentilation(t)
            % If Model 1 re-solved from its own target_VA it would ignore the
            % supplied ventilation; changing target_VA must have no effect.
            c1 = t.cfgM1; c1.patient.target_VA_L_min = 99;
            a = co2.coupling(t.cfgCO2, t.cfgM1, t.cfgM2, 'ETT_7_5');
            b = co2.coupling(t.cfgCO2, c1,      t.cfgM2, 'ETT_7_5');
            t.verifyEqual(b.effort.V_T, a.effort.V_T, 'RelTol', 1e-12);
        end

        function suppliedVentilationIsValidated(t)
            dev = wob.getDevice(t.cfgM1, 'ETT_7_5');
            bad = struct('V_T', 0.4, 'RR', 18);   % missing V_D, V_A, V_E
            t.verifyError(@() wob.simulateModeB(t.cfgM1, dev, t.cfgM1.pattern, Ventilation=bad), ...
                'wob:simulateModeB:incompleteVentilation');
        end

        function couplingProducesALoadModel2CanUse(t)
            cp = co2.coupling(t.cfgCO2, t.cfgM1, t.cfgM2, 'ETT_7_5');
            t.verifyGreaterThan(cp.L_total, 0);
            t.verifyEqual(cp.u_crit_pairing, 0.40);
        end

        function trachLoadRemainsBelowETTLoadThroughTheChain(t)
            LE = co2.coupling(t.cfgCO2, t.cfgM1, t.cfgM2, 'ETT_7_5').L_total;
            LT = co2.coupling(t.cfgCO2, t.cfgM1, t.cfgM2, 'TRACH_8_0').L_total;
            t.verifyLessThan(LT, LE);
        end

        function impliedTargetVAMatchesModel1sConfiguredValue(t)
            % Model 1's target_VA is now COMPUTED here rather than estimated.
            % This asserts the two do not drift apart: if Model 1's config is
            % edited by hand, or this model's constant changes, one of them is
            % stale and the coupling is quietly inconsistent.
            cp = co2.coupling(t.cfgCO2, t.cfgM1, t.cfgM2, 'ETT_7_5', VCO2=300);
            t.verifyEqual(cp.target_VA_implied, t.cfgM1.patient.target_VA_L_min, ...
                'RelTol', 1e-3, ...
                'Model 1 config target_VA has drifted from what Model 1b computes');
        end

        function diseaseGridMatchesTheMetabolicAxis(t)
            % The three severity levels must be the computed V_A at
            % VCO2 = 200/300/400, not hand-scaled approximations.
            expected = arrayfun(@(v) ...
                co2.coupling(t.cfgCO2, t.cfgM1, t.cfgM2, 'ETT_7_5', VCO2=v).target_VA_implied, ...
                [200 300 400]);
            t.verifyEqual(t.cfgM1.disease_grid.target_VA_L_min(:)', expected, 'RelTol', 2e-3);
        end
    end
end
