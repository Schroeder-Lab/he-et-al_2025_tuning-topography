function data = getCircleInfo(folder)

if ~isfile(fullfile(folder, 'circles.times.npy'))
    data = [];
    return
end
t = readNPY(fullfile(folder, 'circles.times.npy'));
data.times = t(:,1);
data.diameter = readNPY(fullfile(folder, 'circles.diameters.npy'));
data.isWhite = readNPY(fullfile(folder, 'circles.isWhite.npy'));
data.xPos = readNPY(fullfile(folder, 'circles.xPos.npy'));
data.yPos = readNPY(fullfile(folder, 'circles.yPos.npy'));
data.interval = readNPY(fullfile(folder, '_ss_recordings.circles_intervals.npy'));