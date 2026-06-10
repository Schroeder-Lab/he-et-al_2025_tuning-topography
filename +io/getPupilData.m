function data = getPupilData(folder)
%GETPUPILDATA   Load pupil data.

% INPUTS
% folder    path to data of recording session

% OUTPUTS
% data
%   .xyPos  [time x 1], coordinates of pupil centre (from top left of video)
%   .time   [time x 1], sampling times of pupil trace

data.xyPos = readNPY(fullfile(folder, 'eye.xyPos.npy'));
data.time = readNPY(fullfile(folder, 'eye.timestamps.npy'));