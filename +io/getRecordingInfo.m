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
%               width), size of reference frame for ROI masks
%   .fovMicrons [nPlanes x 2], size of field-of-view in microns (height x
%               width)
%   .fovBoundaries [nPlanes x 4], for each plane: [top bottom left right]
%               pixels relative to full imaged FOV, pixels outside these
%               boundaries were disregarded for ROI detection as the edges
%               were not always visible
%   .meanFrame  [nPlanes x rows x columns], mean frame after frame
%               alignment
%   .surface    [1], depth of brain surface so that surface - ROI depth =
%               absolute depth for each ROI

data.roiMasks = readNPY(fullfile(folder, '_ss_2pRois.masks.npy'));
data.roiPositions = readNPY(fullfile(folder, '_ss_2pRois.xyz.npy'));
data.fovPix = readNPY(fullfile(folder, '_ss_2pPlanes.fovSizePix.npy'));
data.fovMicrons = readNPY(fullfile(folder, '_ss_2pPlanes.fovSizeMicrons.npy'));
data.fovBoundaries = readNPY(fullfile(folder, '_ss_2pPlanes.fovBoundariesPix.npy'));
data.meanFrame = readNPY(fullfile(folder, '_ss_2pPlanes.meanFrame.npy'));

if isfile(fullfile(folder, '_ss_recordings.brainSurface.npy'))
    data.surface = readNPY(fullfile(folder, '_ss_recordings.brainSurface.npy'));
else
    data.surface = NaN;
end