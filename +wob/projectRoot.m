function root = projectRoot()
%PROJECTROOT Absolute path to the repository root.
%
%   Resolved from this file's own location so that scripts work regardless
%   of the current working folder.

root = fileparts(fileparts(mfilename('fullpath')));
end
