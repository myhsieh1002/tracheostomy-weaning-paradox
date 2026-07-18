classdef tModel2 < matlab.unittest.TestCase
%TMODEL2 Model 2: load-capacity dynamics, folds, coupling, rescue API.
%
%   Covers the sanity checks in section 7 of the Model 2 build plan, plus
%   the two structural claims the analysis rests on -- the scale invariance
%   and the closed-form branch -- each tested against an independent route
%   rather than against itself.

    properties
        cfg1
        cfg2
        p
    end

    methods (TestClassSetup)
        function setup(t)
            addpath(fileparts(fileparts(mfilename('fullpath'))));
            t.cfg1 = wob.loadConfig();
            t.cfg2 = lc.loadConfig();
            t.p    = lc.calibrate(t.cfg2);
        end
    end

    methods (Test)

        % ---------------- spec section 7 limits ----------------

        function noFatigueGivesSingleStablePointAtCmax(t)
            % beta -> 0: recovery alone, so C -> C_max and there is no
            % bistability at any load.
            q = t.p; q.beta = 1e-9;
            C_max = 50;
            fp = lc.fixedPoints(20, C_max, q);
            t.verifyNumElements(fp.C, 1);
            t.verifyEqual(fp.C(1), C_max, 'RelTol', 1e-6);
            t.verifyTrue(fp.stable(1));
            t.verifyFalse(lc.normalizedFolds(q).isBistable);
        end

        function hardThresholdDrivesFoldToUcrit(t)
            % s -> inf: the fold should approach u_crit from below.
            l = zeros(1,4); ss = [100 500 2000 10000];
            for k = 1:numel(ss)
                q = t.p; q.s = ss(k);
                l(k) = lc.normalizedFolds(q).l_high;
            end
            t.verifyTrue(all(diff(l) > 0), 'l_high must rise monotonically with s');
            t.verifyTrue(all(l < t.p.u_crit), 'l_high must stay below u_crit');
            t.verifyEqual(l(end), t.p.u_crit, 'AbsTol', 0.005);
        end

        function middleBranchIsUnstableAndOuterBranchesStable(t)
            C_max = 50;
            nf = lc.normalizedFolds(t.p);
            L = mean([nf.l_low nf.l_high]) * C_max;
            fp = lc.fixedPoints(L, C_max, t.p);
            t.assertNumElements(fp.C, 3);
            t.verifyTrue(fp.stable(1),  'low branch must be stable');
            t.verifyFalse(fp.stable(2), 'middle branch must be unstable');
            t.verifyTrue(fp.stable(3),  'high branch must be stable');
            t.verifyGreaterThan(fp.lambda(2), 0);
        end

        function threeFixedPointsInsideWindowOneOutside(t)
            C_max = 50;
            nf = lc.normalizedFolds(t.p);
            inside = mean([nf.l_low nf.l_high]) * C_max;
            below  = nf.l_low  * C_max * 0.5;
            above  = nf.l_high * C_max * 1.5;
            t.verifyNumElements(lc.fixedPoints(inside, C_max, t.p).C, 3);
            t.verifyNumElements(lc.fixedPoints(below,  C_max, t.p).C, 1);
            t.verifyNumElements(lc.fixedPoints(above,  C_max, t.p).C, 1);
        end

        function capacityIsNeverZeroAtEquilibrium(t)
            % As C -> 0, dC/dt -> alpha*C_max > 0, so C = 0 is not an
            % attractor: the failure state is low-but-positive capacity.
            C_max = 50;
            t.verifyGreaterThan(lc.dCdt(1e-6, 100, C_max, t.p), 0);
            fp = lc.fixedPoints(nan2(lc.normalizedFolds(t.p).l_high)*C_max*2, C_max, t.p);
            t.verifyGreaterThan(min(fp.C), 0);
        end

        function lowBranchDepthMatchesClosedForm(t)
            % C_low/C_max -> alpha/(alpha+beta), which is what pins beta to
            % Laghi's measured capacity drop.
            C_max = 50;
            nf = lc.normalizedFolds(t.p);
            fp = lc.fixedPoints(nf.l_high * C_max * 3, C_max, t.p);
            t.verifyEqual(min(fp.C)/C_max, t.p.alpha/(t.p.alpha+t.p.beta), 'RelTol', 0.02);
        end

        % ---------------- structural claims ----------------

        function systemIsScaleInvariantInCapacity(t)
            % The load-capacity dynamics must depend only on x = C/C_max and
            % l = L/C_max. Verified over a 200-fold range of C_max.
            ref = lc.normalizedFolds(t.p);
            for C_max = [5 20 50 200 1000]
                br = lc.equilibriumBranch(C_max, t.p);
                t.verifyEqual(br.L_fold_low/C_max,  ref.l_low,  'RelTol', 1e-8);
                t.verifyEqual(br.L_fold_high/C_max, ref.l_high, 'RelTol', 1e-8);
            end
        end

        function foldDependsOnlyOnBetaOverAlpha(t)
            % Absolute rates set the timescale, not the bifurcation. This is
            % what retires the "alpha cannot be calibrated" limitation for
            % H1-H3.
            ref = lc.normalizedFolds(t.p);
            ratio = t.p.beta / t.p.alpha;
            for alpha = [0.01 0.05 0.5 1.0]
                q = t.p; q.alpha = alpha; q.beta = ratio*alpha;
                nf = lc.normalizedFolds(q);
                t.verifyEqual(nf.l_low,  ref.l_low,  'RelTol', 1e-8);
                t.verifyEqual(nf.l_high, ref.l_high, 'RelTol', 1e-8);
            end
        end

        function fastFoldsAgreeWithFullBranch(t)
            % lc.normalizedFolds solves the folds directly for speed; it must
            % agree with reading them off the 20k-point branch.
            for s = [30 82.5 200]
                q = t.p; q.s = s;
                nf = lc.normalizedFolds(q);
                br = lc.equilibriumBranch(1, q, PhysicalOnly=false);
                t.verifyEqual(nf.l_low,  br.L_fold_low,  'AbsTol', 1e-9);
                t.verifyEqual(nf.l_high, br.L_fold_high, 'AbsTol', 1e-9);
            end
        end

        function closedFormBranchAgreesWithForwardRootFinding(t)
            % The branch inverts L(C); fixedPoints solves the forward problem
            % by bracketing. They are independent routes to the same curve.
            C_max = 50;
            nf = lc.normalizedFolds(t.p);
            for L = linspace(0.3, nf.l_high*C_max*1.3, 12)
                fp = lc.fixedPoints(L, C_max, t.p);
                for i = 1:numel(fp.C)
                    t.verifyEqual(lc.dCdt(fp.C(i), L, C_max, t.p), 0, 'AbsTol', 1e-8);
                end
            end
        end

        function analyticJacobianMatchesFiniteDifference(t)
            C_max = 50; L = 15;
            for C = [8 20 45]
                [~, dAnalytic] = lc.dCdt(C, L, C_max, t.p);
                h = 1e-6;
                dNumeric = (lc.dCdt(C+h, L, C_max, t.p) - lc.dCdt(C-h, L, C_max, t.p)) / (2*h);
                t.verifyEqual(dAnalytic, dNumeric, 'RelTol', 1e-5);
            end
        end

        function integrationConvergesToPredictedAttractor(t)
            % Spec section 7: a long ode45 integration must land on the
            % attractor the fixed-point analysis predicts.
            %
            % ode45's DEFAULT RelTol is 1e-3, so asserting agreement to 1e-3
            % would be measuring the integrator's own tolerance rather than
            % the model. Tighten the solver well past the assertion.
            C_max = 50;
            nf = lc.normalizedFolds(t.p);
            L = mean([nf.l_low nf.l_high]) * C_max;
            fp = lc.fixedPoints(L, C_max, t.p);
            separatrix = fp.C(2);

            odeOpts = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);
            f = @(tt,C) lc.dCdt(max(C,1e-9), L, C_max, t.p);

            [~, yHi] = ode45(f, [0 800], separatrix*1.05, odeOpts);
            t.verifyEqual(yHi(end), fp.C(3), 'RelTol', 1e-6);

            [~, yLo] = ode45(f, [0 800], separatrix*0.95, odeOpts);
            t.verifyEqual(yLo(end), fp.C(1), 'RelTol', 1e-6);
        end

        function fatigueRateIsStableForExtremeArguments(t)
            % The logistic is evaluated in a branch-free stable form; the
            % naive 1/(1+exp(-x)) overflows here.
            q = t.p; q.s = 500;
            g = lc.fatigueRate([1e-8 0.001 0.4 10 1e6], q);
            t.verifyTrue(all(isfinite(g)));
            t.verifyTrue(all(g >= 0 & g <= 1));
            t.verifyEqual(lc.fatigueRate(q.u_crit, q), 0.5, 'AbsTol', 1e-12);
        end

        % ---------------- validation against published data ----------------

        function reproducesVassilakopoulosDiscrimination(t)
            % The calibration exists to reproduce this; assert it stays
            % reproduced. Pi/Pimax 0.31 weaned, 0.46 failed.
            [~, info] = lc.calibrate(t.cfg2);
            t.verifyTrue(info.discriminates, ...
                'The fold must separate the observed success and failure groups');
            t.verifyTrue(info.predictSuccessWeans);
            t.verifyTrue(info.predictFailureFails);
        end

        function calibrationHitsItsLaghiConstraint(t)
            [p2, info] = lc.calibrate(t.cfg2);
            t.verifyEqual(p2.alpha/(p2.alpha+p2.beta), info.capacityDrop, 'RelTol', 1e-9);
        end

        function unreachableFoldTargetIsRejected(t)
            % l_high is bounded above by u_crit; asking for more must fail
            % loudly rather than silently return a wrong s.
            t.verifyError(@() lc.calibrate(t.cfg2, TargetLHigh=0.45), 'lc:calibrate:noBracket');
        end

        % ---------------- coupling ----------------

        function couplingTakesLoadFromModel1(t)
            % End-to-end: L must come from Model 1, never hard-coded.
            cp = lc.coupling(t.cfg1, t.cfg2, 'ETT_7_5');
            e  = wob.simulateEffort(t.cfg1, 'ETT_7_5');
            t.verifyEqual(cp.L_total, e.P_mus_mean, 'RelTol', 1e-12);
        end

        function trachLoadIsBelowEttLoad(t)
            LE = lc.coupling(t.cfg1, t.cfg2, 'ETT_7_5').L_total;
            LT = lc.coupling(t.cfg1, t.cfg2, 'TRACH_8_0').L_total;
            t.verifyLessThan(LT, LE);
        end

        function utilisationPairsThresholdWithConvention(t)
            % Mixing a duty-weighted load with the mean-pressure threshold
            % would mis-scale the whole bifurcation; the pairing is enforced.
            e = wob.simulateEffort(t.cfg1, 'ETT_7_5');
            [Lm, mm] = lc.utilisation(e, 'mean');
            [Ld, md] = lc.utilisation(e, 'duty_weighted');
            t.verifyEqual(mm.u_crit_pairing, 0.40);
            t.verifyEqual(md.u_crit_pairing, 0.15);
            t.verifyEqual(Ld, Lm * e.dutyCycle, 'RelTol', 1e-9);

            % 'peak' has no published threshold, so it must refuse to supply
            % one rather than quietly hand back a plausible-looking number.
            [Lp, mp] = lc.utilisation(e, 'peak');
            t.verifyEqual(Lp, e.P_mus_peak, 'RelTol', 1e-12);
            t.verifyTrue(isnan(mp.u_crit_pairing));

            t.verifyError(@() lc.utilisation(e,'rms'), 'lc:utilisation:unknownConvention');
        end

        function couplingApportionmentIsOptional(t)
            cheap = lc.coupling(t.cfg1, t.cfg2, 'ETT_7_5');
            full  = lc.coupling(t.cfg1, t.cfg2, 'ETT_7_5', struct(), Apportionment=true);
            t.verifyTrue(isnan(cheap.L_device_direct));
            t.verifyFalse(isnan(full.L_device_direct));
            t.verifyEqual(cheap.L_total, full.L_total, 'RelTol', 1e-12);
        end

        % ---------------- rescue API (spec section 2.5b) ----------------

        function rescueWindowReturnsSpecifiedFields(t)
            rw = lc.rescueWindow(50, [5 8], t.p);
            for f = ["fold_left","fold_right","bistable","window","basin_at","C_max_window","spans_fold"]
                t.verifyTrue(isfield(rw, f), sprintf('rescueWindow must expose %s', f));
            end
            t.verifyClass(rw.basin_at, 'function_handle');
        end

        function rescueWindowFoldsScaleWithCapacity(t)
            a = lc.rescueWindow(30, [5 8], t.p);
            b = lc.rescueWindow(60, [5 8], t.p);
            t.verifyEqual(b.fold_right, 2*a.fold_right, 'RelTol', 1e-9);
        end

        function rescueWindowAgreesWithRescueOutcome(t)
            % Two implementations of the same rule must not drift apart.
            nf = lc.normalizedFolds(t.p);
            for C_max = [20 35 50 70]
                LE = nf.l_high*C_max*1.1; LT = nf.l_high*C_max*0.9;
                r  = lc.rescueOutcome(LE, LT, C_max, t.p);
                rw = lc.rescueWindow(C_max, [LT LE], t.p);
                t.verifyEqual(r.isRescued, rw.spans_fold);
                t.verifyEqual(r.L_fold, rw.fold_right, 'RelTol', 1e-12);
            end
        end

        function rescueOutcomeClassifiesAllThreeRegimes(t)
            C_max = 50;
            nf = lc.normalizedFolds(t.p);
            Lf = nf.l_high * C_max;
            t.verifyEqual(lc.rescueOutcome(Lf*0.8, Lf*0.7, C_max, t.p).outcome, "both_wean");
            t.verifyEqual(lc.rescueOutcome(Lf*1.1, Lf*0.9, C_max, t.p).outcome, "rescued");
            t.verifyEqual(lc.rescueOutcome(Lf*1.5, Lf*1.3, C_max, t.p).outcome, "both_fail");
        end

        function rescueWindowWidthMatchesClosedForm(t)
            % width = dL_device / l_high
            nf = lc.normalizedFolds(t.p);
            LT = 6; LE = 9;
            rw = lc.rescueWindow(50, [LT LE], t.p);
            t.verifyEqual(diff(rw.C_max_window), (LE-LT)/nf.l_high, 'RelTol', 1e-9);
        end

        function basinAtResolvesTheSeparatrix(t)
            C_max = 50;
            nf = lc.normalizedFolds(t.p);
            L = mean([nf.l_low nf.l_high]) * C_max;
            rw = lc.rescueWindow(C_max, [L L], t.p);
            sep = lc.fixedPoints(L, C_max, t.p).C(2);
            t.verifyEqual(rw.basin_at(sep*1.05, L), "high");
            t.verifyEqual(rw.basin_at(sep*0.95, L), "low");
        end

        function rescueWindowRejectsBadRange(t)
            t.verifyError(@() lc.rescueWindow(50, [9 5], t.p), 'lc:rescueWindow:badRange');
        end

        function trajectoryModeTracksAClosingWindow(t)
            % The integration seam for Model 2b. 2b does not exist yet, so
            % this exercises it against a synthetic decaying capacity.
            tt = linspace(0, 20, 60);
            C_max_t = linspace(70, 15, 60);
            traj = lc.rescueWindowTrajectory(tt, C_max_t, [6 9], t.p);
            t.verifyTrue(all(diff(traj.fold_right) < 0), 'fold must fall as capacity falls');
            t.verifyTrue(traj.windowOpen(1), 'window should start open at high capacity');
            t.verifyFalse(traj.windowOpen(end), 'window should be shut at low capacity');
            t.verifyTrue(isfinite(traj.firstClose));
        end

        function trajectoryRejectsMismatchedLengths(t)
            t.verifyError(@() lc.rescueWindowTrajectory(1:5, 1:4, [6 9], t.p), ...
                'lc:rescueWindowTrajectory:sizeMismatch');
        end

        % ---------------- stochastic / early warning ----------------

        function noiseScalesWithCapacity(t)
            % sigma is a fraction of C_max, so the noise-to-signal ratio is
            % invariant the way the rest of the model is.
            tt = linspace(0, 5, 500);
            s1 = lc.simulateStochastic(tt, repmat(5,1,500), 25, t.p, Sigma=0.02);
            s2 = lc.simulateStochastic(tt, repmat(10,1,500), 50, t.p, Sigma=0.02);
            t.verifyEqual(s2.sigma_abs, 2*s1.sigma_abs, 'RelTol', 1e-12);
        end

        function autocorrelationRisesApproachingTheFold(t)
            % Critical slowing down, measured at a sampling interval
            % comparable to the relaxation time (1/alpha ~ 6.7 h). Sampling
            % at the integration step instead would pin AR(1) at ~0.998 in
            % both arms and the test would compare noise.
            %
            % Noise is kept low enough that neither arm escapes to the low
            % branch: once collapsed, the state sits where the Jacobian is
            % steeply negative, relaxation is FAST, and AR(1) drops -- which
            % would invert this comparison for a reason that has nothing to
            % do with critical slowing down. The no-collapse precondition is
            % asserted rather than assumed.
            C_max = 50;
            nf = lc.normalizedFolds(t.p);
            tt = linspace(0, 800, 160000);
            subN = round(2.0/(tt(2)-tt(1)));
            sigma = 0.004;

            far  = repmat(nf.l_low*C_max*1.05,  1, numel(tt));
            near = repmat(nf.l_high*C_max*0.99, 1, numel(tt));

            sFar  = lc.simulateStochastic(tt, far,  C_max, t.p, Sigma=sigma, Seed=7);
            sNear = lc.simulateStochastic(tt, near, C_max, t.p, Sigma=sigma, Seed=7);

            sepNear = lc.fixedPoints(near(1), C_max, t.p).C(2);
            t.assumeTrue(all(sNear.C > sepNear), ...
                'precondition: the near-fold path must not collapse, or AR(1) reports the low branch');

            aFar  = mean(lc.earlyWarning(tt, sFar.C,  Subsample=subN, Window=100).ar1, 'omitnan');
            aNear = mean(lc.earlyWarning(tt, sNear.C, Subsample=subN, Window=100).ar1, 'omitnan');
            t.verifyGreaterThan(aNear, aFar);
        end

        function stochasticRejectsNonMonotonicTime(t)
            t.verifyError(@() lc.simulateStochastic([0 2 1], [5 5 5], 50, t.p), ...
                'lc:simulateStochastic:nonMonotonicTime');
        end

        % ---------------- config integrity ----------------

        function unknownDynamicalParameterIsRejected(t)
            t.verifyError(@() lc.params(t.cfg2, struct('gamma',1)), 'lc:params:unknownParameter');
        end

        function capacityGridCoversTheObservedWeaningRange(t)
            % Vassilakopoulos: MIP 42.3 (failure) / 53.8 (success). The grid
            % must bracket the band where failure actually happens.
            g = t.cfg2.capacity_grid.C_max;
            t.verifyLessThanOrEqual(min(g), 30, 'grid must reach the weak end where failure lives');
            t.verifyGreaterThanOrEqual(max(g), 60);
        end
    end
end

function y = nan2(x)
y = x; if isnan(y), y = 1; end
end
