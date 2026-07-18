function runAll(opts)
%RUNALL One-command reproduction of every result in the repository.
%
%   runAll()                 run tests, all experiments, regenerate summaries
%   runAll(Tests=false)      skip the test suite
%   runAll(Only="m2")        one group: "m1" | "m1b" | "m2" | "m2b" | "all"
%   runAll(Visible=true)     show figures as they are built
%
%   Headless (does not touch a shared MATLAB desktop):
%
%       matlab -batch "runAll"
%
%   Outputs land in results/figures (PNG at 300 dpi + vector PDF),
%   results/tables (CSV), and results/summary.md + results/summary_m2.md.

arguments
    opts.Tests   (1,1) logical = true
    opts.Only    (1,1) string  = "all"
    opts.Visible (1,1) logical = false
end

root = fileparts(mfilename('fullpath'));
addpath(root);
addpath(fullfile(root, 'experiments'));
addpath(fullfile(root, 'manuscript'));

for d = ["results/figures", "results/tables"]
    p = fullfile(root, d);
    if ~isfolder(p), mkdir(p); end
end

prevVis = get(0, 'DefaultFigureVisible');
set(0, 'DefaultFigureVisible', matlab.lang.OnOffSwitchState(opts.Visible));
warning('off', 'MATLAB:print:ContentTypeImageSuggested');
cleanup = onCleanup(@() set(0, 'DefaultFigureVisible', prevVis));

t0 = tic;

if opts.Tests
    banner('TESTS');
    r = [runtests(fullfile(root,'tests','tModel1.m')), ...
         runtests(fullfile(root,'tests','tModel2.m')), ...
         runtests(fullfile(root,'tests','tModel1b.m')), ...
         runtests(fullfile(root,'tests','tModel2b.m'))];
    fprintf('\n  %d passed, %d failed, %d incomplete\n', ...
        sum([r.Passed]), sum([r.Failed]), sum([r.Incomplete]));
    if any([r.Failed])
        error('runAll:testsFailed', ...
            ['%d test(s) failed. Results are not regenerated from a failing ' ...
             'model -- fix the failures first.'], sum([r.Failed]));
    end
end

R = struct();

if ismember(opts.Only, ["all","m1"])
    banner('MODEL 1');
    step('exp1  resistance curves (H1)');        R.exp1 = exp1_resistanceCurves();
    step('exp2  imposed WOB / PTP (H2)');        R.exp2 = exp2_imposedWOB();
    step('exp3  dead space and effort (H3)');    R.exp3 = exp3_deadspaceEffort();
    step('exp4  disease severity grid (H4)');    R.exp4 = exp4_diseaseSeverityGrid();
    step('exp5  global sensitivity');            R.exp5 = exp5_sensitivity();
end

if ismember(opts.Only, ["all","m1b"])
    banner('MODEL 1b -- CO2 KINETICS');
    step('co2_exp1  steady-state PaCO2 (H1)');       R.co2_exp1 = co2_exp1_steadyState();
    step('co2_exp2  required V_E, and its cost (H2)'); R.co2_exp2 = co2_exp2_requiredVE();
    step('co2_exp3  dilution grid (H3)');            R.co2_exp3 = co2_exp3_dilutionGrid();
    step('co2_exp4  transient after tracheostomy');  R.co2_exp4 = co2_exp4_transient();
end

if ismember(opts.Only, ["all","m2"])
    banner('MODEL 2');
    step('m2_exp1  bifurcation + hysteresis (H1)');  R.m2_exp1 = m2_exp1_bifurcation1p();
    step('m2_exp2  two-parameter structure (H2)');   R.m2_exp2 = m2_exp2_cusp2p();
    step('m2_exp3  rescue window (H3)');             R.m2_exp3 = m2_exp3_rescueWindow();
    step('m2_exp4  early warning (H4)');             R.m2_exp4 = m2_exp4_earlyWarning();
    step('m2_exp5  sensitivity + invariance test');  R.m2_exp5 = m2_exp5_sensitivity();
end

if ismember(opts.Only, ["all","m2b"])
    banner('MODEL 2b -- VIDD CAPACITY DYNAMICS');
    step('vidd_exp1  capacity vs activity, both shapes (H1)'); R.vidd_exp1 = vidd_exp1_shape();
    step('vidd_exp2  capacity trajectories (H2)');             R.vidd_exp2 = vidd_exp2_trajectory();
    step('vidd_exp3  catabolic: offset vs drain (H3)');        R.vidd_exp3 = vidd_exp3_catabolic();
    step('vidd_exp4  dynamic rescue window (H4, capstone)');   R.vidd_exp4 = vidd_exp4_dynamicRescue();
    step('vidd_exp5  the A-gap, stress-tested');               R.vidd_exp5 = vidd_exp5_Agap();
end

if ismember(opts.Only, "all")
    banner('MANUSCRIPT');
    step('Figure 1  four-model schematic');   make_fig1();
end

banner('SUMMARIES');
writeSummary(R, opts.Only);

fprintf('\nDone in %.0f s. See results/summary.md and results/summary_m2.md\n', toc(t0));
end

function banner(s)
fprintf('\n%s\n=== %s\n%s\n', repmat('=',1,64), s, repmat('=',1,64));
end

function step(s)
fprintf('\n--- %s\n', s);
end
