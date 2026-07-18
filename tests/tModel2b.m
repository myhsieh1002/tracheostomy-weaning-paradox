classdef tModel2b < matlab.unittest.TestCase
%TMODEL2B Model 2b: VIDD capacity dynamics and the dynamic coupling.
%
%   Covers section 6 of the Model 2b build plan, the structural claims the
%   model rests on (monostability, the shape being g's, the slow-fast
%   separation), and the calibration against Jaber / Zambon / Lecronier.

    properties
        cfgV
        cfgM1
        cfgM2
        p
        info
    end

    methods (TestClassSetup)
        function setup(t)
            addpath(fileparts(fileparts(mfilename('fullpath'))));
            t.cfgV  = vidd.loadConfig();
            t.cfgM1 = wob.loadConfig();
            t.cfgM2 = lc.loadConfig();
            [t.p, t.info] = vidd.calibrate(t.cfgV);
        end
    end

    methods (Test)

        % ---------------- structure ----------------

        function systemIsMonostable(t)
            % The defining difference from Model 2: A is exogenous, the
            % system is linear in C, so there is exactly one equilibrium and
            % no fold, at any parameters. This is what makes the coupling a
            % clean slow-fast one.
            for A = [0 0.3 0.6 1]
                [~, dfdC] = vidd.dCdt(50, A, t.p);
                t.verifyLessThan(dfdC, 0, sprintf('equilibrium must be stable at A=%g', A));
            end
        end

        function equilibriumMatchesLongIntegration(t)
            % The closed-form C* must equal where a long integration lands.
            A = 0.4;
            Cstar = vidd.equilibriumCapacity(A, t.p);
            r = vidd.simulateCapacity(t.cfgV, t.p, @(tt) A, C0=10, Days=200);
            t.verifyEqual(r.C(end), Cstar, 'RelTol', 1e-3);
        end

        function analyticJacobianMatchesFiniteDifference(t)
            A = 0.4;
            for C = [10 40 70]
                [~, dAn] = vidd.dCdt(C, A, t.p);
                h = 1e-6;
                dNum = (vidd.dCdt(C+h,A,t.p) - vidd.dCdt(C-h,A,t.p))/(2*h);
                t.verifyEqual(dAn, dNum, 'RelTol', 1e-6);
            end
        end

        function capacityFactorEntersDynamics(t)
            % The (C_max0 - C) factor the spec omitted: C_max0 must affect
            % the trajectory, not just the initial condition. Doubling it
            % must change C* (it would not if the term were k_syn*h alone).
            a = t.p; a.C_max0 = 80;
            b = t.p; b.C_max0 = 160;
            t.verifyNotEqual(vidd.equilibriumCapacity(0.5, a), vidd.equilibriumCapacity(0.5, b));
        end

        function equilibriumIsAFractionOfCeiling(t)
            % C* = syn/(syn+loss) * C_max0, so C* < C_max0 always and the
            % fraction is scale-free in C_max0.
            for A = [0 0.5 1]
                frac = vidd.equilibriumCapacity(A, t.p) / t.p.C_max0;
                t.verifyGreaterThanOrEqual(frac, 0);
                t.verifyLessThanOrEqual(frac, 1 + 1e-9);
            end
        end

        % ---------------- the shape is g's, not h's ----------------

        function shapeSurvivesConstantSynthesis(t)
            % Set h0=1 (h constant); the shape of C*(A) must be unchanged in
            % CHARACTER -- monotone stays monotone, U stays U -- proving the
            % shape is the degradation term's.
            A = linspace(0,1,200);
            for mode = ["monotonic","ushape"]
                q1 = t.p; q1.g_mode = mode;
                q2 = q1;  q2.h0 = 1;
                c1 = vidd.equilibriumCapacity(A, q1);
                c2 = vidd.equilibriumCapacity(A, q2);
                % number of sign changes in the derivative = shape signature
                t.verifyEqual(localTurns(c1), localTurns(c2), ...
                    sprintf('%s: constant h must not change the shape', mode));
            end
        end

        function monotonicModeIsMonotone(t)
            A = linspace(0,1,200);
            C = vidd.equilibriumCapacity(A, setMode(t.p,"monotonic"));
            t.verifyTrue(all(diff(C) > 0), 'monotonic mode: C* must rise with A');
        end

        function ushapeModeHasAnInteriorPeak(t)
            A = linspace(0,1,400);
            q = setMode(t.p,"ushape");
            C = vidd.equilibriumCapacity(A, q);
            [~, iMax] = max(C);
            t.verifyGreaterThan(iMax, 1);
            t.verifyLessThan(iMax, numel(A));
        end

        function degradationFloorIsOneAtZeroActivity(t)
            % Convention g(0)=1 so k_deg is the complete-disuse rate.
            t.verifyEqual(vidd.degradation(0, setMode(t.p,"monotonic")), 1, 'AbsTol', 1e-9);
        end

        function unknownGModeErrors(t)
            q = t.p; q.g_mode = "sinusoidal";
            t.verifyError(@() vidd.degradation(0.5, q), 'vidd:degradation:unknownMode');
        end

        % ---------------- calibration ----------------

        function kDegIsAnchoredToJaber(t)
            t.verifyEqual(t.info.k_deg, 0.064, 'AbsTol', 1e-6);
        end

        function reproducesJaberForceLoss(t)
            % Independent check: complete disuse for 6 days should lose ~32%,
            % which is Jaber's measured force loss.
            t.verifyEqual(t.info.jaber_predicted_loss_at_A0, 0.32, 'AbsTol', 0.03);
        end

        function kSynIsIdentifiedFromLecronier(t)
            % The two-group solve must succeed and yield a positive rate.
            t.verifyTrue(t.info.k_syn_identified);
            t.verifyGreaterThan(t.info.k_syn, 0);
            t.verifyGreaterThan(t.info.lecronier_E, 0);
            t.verifyLessThan(t.info.lecronier_E, 1);
        end

        function lecronierAndDemouleAgreeOnSepsisOffset(t)
            % Two studies, two methods: the offsets must agree to a few
            % percentage points, or the two-group solve is inconsistent with
            % Demoule's regression.
            t.verifyLessThan(t.info.offset_agreement, 0.06, ...
                'Lecronier-derived and Demoule sepsis offsets disagree by >6 pp');
        end

        function zambonFitIsGood(t)
            % Normalised by Zambon's own CMV rate, the monotone form should
            % fit the three negative points well.
            t.verifyLessThan(t.info.zambon_rmse, 0.12);
        end

        function scaleConflictIsReportedNotHidden(t)
            % Jaber and Zambon disagree ~22% on the disuse rate; the code
            % must surface that rather than bury it in the fit.
            t.verifyGreaterThan(t.info.scale_conflict_ratio, 1.1);
            t.verifyLessThan(t.info.scale_conflict_ratio, 1.4);
        end

        % ---------------- sepsis: offset vs rate ----------------

        function offsetModeLowersInitialCapacityNotRate(t)
            % d_mode=offset: septic and non-septic share dynamics, differ
            % only in C(0).
            [c0n, ~] = vidd.initialCapacity(t.p, false);
            [c0s, ~] = vidd.initialCapacity(t.p, true);
            t.verifyEqual(c0n, t.p.C_max0, 'AbsTol', 1e-9);
            t.verifyLessThan(c0s, c0n);
            t.verifyEqual(c0s/c0n, 1 - t.p.sepsis_offset_fraction, 'RelTol', 1e-9);
        end

        function septicPatientRecoversUnderOffsetMode(t)
            % The Lecronier finding the model must reproduce: a septic
            % patient on spontaneous breathing IMPROVES, not declines.
            A = vidd.supportToActivity(t.cfgV.strategy.scenarios.trach_spontaneous.support_level, t.cfgV);
            r = vidd.simulateCapacity(t.cfgV, t.p, @(tt) A, IsSeptic=true, Days=10);
            t.verifyGreaterThan(r.C(end), r.C(1), ...
                'septic patient on spontaneous breathing must recover (Lecronier 2022)');
        end

        function rateModeCollapsesCapacity(t)
            % d_mode=rate (the spec): a real drain drives C down and, at high
            % enough d, no activity holds it.
            q = t.p; q.d_mode = "rate"; q.d_disease = 0.2;
            r = vidd.simulateCapacity(t.cfgV, q, @(tt) 0.6, Days=14);
            t.verifyLessThan(r.C(end), r.C(1));
        end

        % ---------------- strategy leverage (H3) ----------------

        function absoluteLeverageCollapsesWithCatabolism(t)
            % H3, honestly: the ABSOLUTE gap between best and worst strategy
            % shrinks as d_disease rises.
            q0 = t.p; q0.d_mode = "rate"; q0.d_disease = 0.0;
            q1 = t.p; q1.d_mode = "rate"; q1.d_disease = 0.3;
            t.verifyGreaterThan(vidd.strategyLeverage(q0).absolute, ...
                                 vidd.strategyLeverage(q1).absolute);
        end

        % ---------------- coupling to Model 2 ----------------

        function couplingUsesModel2NamedAPI(t)
            % The seam: the trajectory must come back with the fields
            % lc.rescueWindowTrajectory produces, i.e. the fold was solved by
            % Model 2, not re-derived here.
            cp = vidd.couplingToModel2(t.cfgV, t.cfgM1, t.cfgM2, 'trach_spontaneous');
            for f = ["windowOpen","fold_right","spans_fold","firstClose"]
                t.verifyTrue(isfield(cp.traj, f), sprintf('trajectory missing %s', f));
            end
        end

        function couplingTakesLoadFromModel1(t)
            % Loads must come from Model 1, never hard-coded.
            cp = vidd.couplingToModel2(t.cfgV, t.cfgM1, t.cfgM2, 'trach_spontaneous');
            LE = lc.coupling(t.cfgM1, t.cfgM2, 'ETT_7_5').L_total;
            LT = lc.coupling(t.cfgM1, t.cfgM2, 'TRACH_8_0').L_total;
            t.verifyEqual(cp.L_ETT, LE, 'RelTol', 1e-12);
            t.verifyEqual(cp.L_TRACH, LT, 'RelTol', 1e-12);
        end

        function timescaleSeparationIsAdequate(t)
            % The quasi-static coupling needs the capacity axis (days) well
            % slower than Model 2's fatigue axis (hours). Assert it, since
            % the whole coupling rests on it.
            cp = vidd.couplingToModel2(t.cfgV, t.cfgM1, t.cfgM2, 'trach_spontaneous');
            t.verifyGreaterThan(cp.separation, 5, ...
                'timescale separation below 5x: quasi-static coupling not safe');
        end

        function trajectoryFoldFallsAsCapacityFalls(t)
            % A declining capacity must drag Model 2's fold down with it.
            tt = linspace(0, 20, 50);
            Cdrop = linspace(70, 20, 50);
            traj = lc.rescueWindowTrajectory(tt, Cdrop, [6 9], lc.calibrate(t.cfgM2));
            t.verifyTrue(all(diff(traj.fold_right) < 0));
        end

        function unknownScenarioErrors(t)
            t.verifyError(@() vidd.couplingToModel2(t.cfgV, t.cfgM1, t.cfgM2, 'nonesuch'), ...
                'vidd:couplingToModel2:unknownScenario');
        end

        % ---------------- input validation ----------------

        function activityOutOfRangeErrors(t)
            t.verifyError(@() vidd.dCdt(50, 1.5, t.p), 'vidd:dCdt:activityOutOfRange');
            t.verifyError(@() vidd.supportToActivity(1.5, t.cfgV), 'vidd:supportToActivity:outOfRange');
        end

        function unknownParameterOverrideErrors(t)
            t.verifyError(@() vidd.params(t.cfgV, struct('kdeg', 1)), 'vidd:params:unknownParameter');
        end

        function supportMapsToActivityInverted(t)
            % support 1 -> A 0 (fully controlled), support 0 -> A A_max.
            t.verifyEqual(vidd.supportToActivity(1, t.cfgV), 0, 'AbsTol', 1e-12);
            t.verifyEqual(vidd.supportToActivity(0, t.cfgV), t.cfgV.strategy.A_max, 'AbsTol', 1e-12);
        end
    end
end

function n = localTurns(y)
% Number of interior sign changes of the derivative = shape signature.
d = sign(diff(y));
d = d(d ~= 0);
n = sum(abs(diff(d)) > 0);
end

function q = setMode(p, mode)
q = p; q.g_mode = mode;
end
