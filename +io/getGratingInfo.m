function data = getGratingInfo(folder, type)

if nargin < 2
    type = '';
end

switch type
    case 'drifting'
        str = 'gratingsDrifting';
    case 'static'
        str = 'gratingsStatic';
    otherwise
        str = 'grating';
end

if isfile(fullfile(folder, sprintf('_ss_%s.intervals.npy', str)))
    data.times = readNPY(fullfile(folder, ...
        sprintf('_ss_%s.intervals.npy', str)));
else
    data = [];
    return
end
data.ids = readNPY(fullfile(folder, sprintf('_ss_%s._ss_%sID.npy', str, str)));
data.directions = readNPY(fullfile(folder, ...
    sprintf('_ss_%sID.directions.npy', str)));
data.interval = readNPY(fullfile(folder, ...
    sprintf('_ss_recordings.%s_intervals.npy', str)));