function differences = determinePreferenceDiff(angles, type)
%DETERMINEPREFERENCEDIFF   Calculate all pairwise differences in preferred
%directions/orientations.

% INPUTS
% angles        [ROIs x 1], preferred directions or orientations of ROIs
% type          'ori' or 'dir'

% OUTPUTS
% differences   [pairs x 1], pairwise differences

differences = [];
diff = abs(angles - angles');
diff = diff(tril(true(size(diff)),-1));
if strcmp(type, 'dir')
    ind = diff > 180;
    diff(ind) = 360 - diff(ind);
elseif  strcmp(type, 'ori')
    ind = diff > 90;
    diff(ind) = 180 - diff(ind);
else
    return
end
differences = diff;