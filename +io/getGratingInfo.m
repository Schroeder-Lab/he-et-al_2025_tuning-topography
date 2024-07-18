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
data.interval = readNPY(fullfile(folder, ...
    sprintf('_ss_recordings.%s_intervals.npy', type)));
switch type
    case {'gratingsDrifting', 'bars'}
        data.directions = readNPY(fullfile(folder, ...
            sprintf('_ss_%sID.directions.npy', type)));
    case 'gratingsStatic'
        data.orientations = readNPY(fullfile(folder, ...
            sprintf('_ss_%sID.orientations.npy', type)));
        data.phases = readNPY(fullfile(folder, ...
            sprintf('_ss_%sID.spatialPhases.npy', type)));
    otherwise
        data.directions = readNPY(fullfile(folder, ...
            sprintf('_ss_%sID.directions.npy', type)));
end