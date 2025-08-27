function data = getVisNoiseInfo(folder)
%GETVISNOISEINFO   Load visual noise stimulus data.

% INPUTS
% folder        path to data of recording session

% OUTPUTS
% data
%   .stimOrder  [frameIDs], IDs of stimulus frame presented in that order
%   .times      [frameIDs], time of frame presentations
%   .edges      [left right top bottom], edges of stimulus frames in visual
%               degree
%   .frames     [frames x rows x columns], pixel values of each unique
%               stimulus frame

if ~isfile(fullfile(folder, '_ss_sparseNoise.times.npy'))
    data = [];
    return
end
data.stimOrder = readNPY(fullfile(folder, '_ss_sparseNoise._ss_sparseNoiseID.npy'));
data.times = readNPY(fullfile(folder, '_ss_sparseNoise.times.npy'));
data.edges = readNPY(fullfile(folder, '_ss_sparseNoiseArea.edges.npy'));
data.frames = readNPY(fullfile(folder, '_ss_sparseNoiseID.map.npy'));