function label = deviceLabel(key)
%DEVICELABEL Human-readable axis label for a device config key.
%
%   label = viz.deviceLabel('ETT_7_5') returns 'ETT 7.5 mm'.
%
%   Config keys encode the internal diameter with underscores because a
%   struct field cannot contain a decimal point. A plain strrep of '_' to
%   ' ' therefore renders 'ETT 7 5', which is wrong in a published figure.
%   Keys with no numeric size ('NATIVE_UPPER_AIRWAY') fall back to
%   lower-case words.

arguments
    key (1,:) char
end

tok = strsplit(key, '_');
isNum = ~cellfun(@isempty, regexp(tok(2:end), '^\d+$', 'once'));

if numel(tok) > 1 && all(isNum)
    label = sprintf('%s %s mm', tok{1}, strjoin(tok(2:end), '.'));
else
    label = lower(strrep(key, '_', ' '));
end
end
