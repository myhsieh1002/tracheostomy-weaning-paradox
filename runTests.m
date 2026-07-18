function runTests(which)
%RUNTESTS Run the test suites and print a compact summary.
%
%   runTests()          run every suite
%   runTests("m1b")     run one suite: "m1" | "m2" | "m1b" | "all"
%
%   Written to be driven headlessly:
%
%       matlab -batch "runTests('m1b')"
%
%   which matters because the MATLAB MCP server attaches to a SHARED desktop
%   instance -- one workspace, one current folder, one function cache across
%   every client. Two Claude Code sessions on it would clobber each other's
%   variables and cwd, and the `clear functions` this repo needs after every
%   edit would wipe the other session's state too. `matlab -batch` starts an
%   independent headless process instead, so runs here cannot disturb a
%   desktop session doing something else. Cost is ~14 s of startup per
%   invocation, which is why this batches.

arguments
    which (1,1) string = "all"
end

root = fileparts(mfilename('fullpath'));
addpath(root);

suites = struct('m1',  'tModel1.m', ...
                'm2',  'tModel2.m', ...
                'm1b', 'tModel1b.m', ...
                'm2b', 'tModel2b.m');

if which == "all"
    names = fieldnames(suites)';
else
    if ~isfield(suites, which)
        error('runTests:unknownSuite', ...
            'Unknown suite "%s". Valid: m1, m2, m1b, all.', which);
    end
    names = {char(which)};
end

r = [];
for k = 1:numel(names)
    f = fullfile(root, 'tests', suites.(names{k}));
    if ~isfile(f)
        fprintf('  (skipping %s -- not present)\n', suites.(names{k}));
        continue;
    end
    r = [r, runtests(f)]; %#ok<AGROW>
end

fprintf('\n===== %d passed, %d failed, %d incomplete of %d =====\n', ...
    sum([r.Passed]), sum([r.Failed]), sum([r.Incomplete]), numel(r));

bad = find([r.Failed] | [r.Incomplete]);
for k = bad
    fprintf('  >> %s\n', r(k).Name);
end

if ~isempty(bad)
    fprintf('\n%d test(s) need attention.\n', numel(bad));
end
end
