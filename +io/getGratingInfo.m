function data = getGratingInfo(folder, type)

if nargin < 2
    type = 'grating';
end

if isfile(fullfile(folder, sprintf('_ss_%s.intervals.npy', type)))
    data.times = readNPY(fullfile(folder, ...
        sprintf('_ss_%s.intervals.npy', type)));
else
    data = [];
    return
end
data.ids = readNPY(fullfile(folder, sprintf('_ss_%s._ss_%sID.npy', type, type)));
data.directions = readNPY(fullfile(folder, ...
    sprintf('_ss_%sID.directions.npy', type)));
data.interval = readNPY(fullfile(folder, ...
    sprintf('_ss_recordings.%s_intervals.npy', type)));