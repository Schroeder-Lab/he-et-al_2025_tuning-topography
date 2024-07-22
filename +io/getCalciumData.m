function data = getCalciumData(folder)
%GETCALCIUMDATA   Load calcium data.

% INPUTS
% folder    path to data of recording session

% OUTPUTS
% data
%   .planes [ROIs x 1], imaging plane of each ROI
%   .ids    [ROIs x 1, ID of each ROI
%   .traces [time x ROIs], calcium trace of each ROI
%   .time   [time x 1], sampling times of calcium traces (plane 1)
%   .delays [plane x 1], temporal delay relative to sampling time (of plane 1)

% for each ROI: its imaging plane, its ID
data.planes = readNPY(fullfile(folder, '_ss_2pRois._ss_2pPlanes.npy'));
data.ids = readNPY(fullfile(folder, '_ss_2pRois.ids.npy'));
% calcium traces of all ROIs and sampling time
data.traces = readNPY(fullfile(folder, '_ss_2pCalcium.dff.npy'));
data.time = readNPY(fullfile(folder, '_ss_2pCalcium.timestamps.npy'));
% for each imaging plane: temporal delay relative to sampling time (actual
% time for first plane)
data.delays = readNPY(fullfile(folder, '_ss_2pPlanes.delay.npy'));