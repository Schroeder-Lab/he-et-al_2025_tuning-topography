function data = getCircleInfo(folder)
%GETCIRCLEINFO   Load circle stimulus data.

% INPUTS
% folder        path to data of recording session

% OUTPUTS
% data
%   .times      [t], times of each stimulus frame presentation in sec
%   .diameter   [t], diameter of presented circle in each stimulus frame
%   .isWhite    [t], gray value of presented circle in each stimulus frame,
%               0: black, 1: white
%   .xPos       [t], azimuth of presented circle in each stimulus frame
%   .yPos       [t], elevation of presented circle in each stimulus frame
%   .interval   [2], start and end of circle stimulus paradigm

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