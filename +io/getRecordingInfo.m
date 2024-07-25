function data = getRecordingInfo(folder)
%GETRECORDINGINFO   Load information about the recording.

% INPUTS
% folder    path to data of recording session

% OUTPUTS
% data
%   .roiMasks   [ROIs x pix], pixel indices with FOV (in pixels) of ROI
%               masks
%   .roiPositions [ROIs x 3], (horizontal, vertical, depth)-coordinates in
%               microns
%   .fovPix     [nPlanes x 2], size of field-of-view in pixels (height x
%               width)
%   .fovMicrons [nPlanes x 2], size of field-of-view in microns (height x
%               width)

data.roiMasks = readNPY(fullfile(folder, '_ss_2pRois.masks.npy'));
data.roiPositions = readNPY(fullfile(folder, '_ss_2pRois.xyz.npy'));
data.fovPix = readNPY(fullfile(folder, '_ss_2pPlanes.fovSizePix.npy'));
data.fovMicrons = readNPY(fullfile(folder, '_ss_2pPlanes.fovSizeMicrons.npy'));