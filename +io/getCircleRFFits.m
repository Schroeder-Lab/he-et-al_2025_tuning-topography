function data = getCircleRFFits(folder)
%GETCIRCLERFFITS   Load receptive fields fitted to responses to circle
%stimuli.

% INPUTS
% folder            path to data of recording session

% OUTPUTS
% data
%   .maps           [ROIs x rows x columns x t_rf x type], spatio-temporal
%                   RF for ON (type=1) and OFF (type=2) subfields
%   .bestSubFields  [ROIs], 1: 'ON', 2: 'OFF', or 3: 'ON+OFF'
%   .subFieldSigns  [ROIs x 2], signs of ON and OFF fields
%   .gaussPars      [ROIs x parameters], (amplitude, xCenter, xSTD,
%                   yCenter, ySTD, rotation, offset)
%   .peaks          [ROIs], amplitude of RF peak (of best type) in STD of
%                   noise
%   .gaussMasks     [ROIs x rows x columns x type], best fitting Gaussian
%                   mask for ON (type=1) and OFF (type=2) subfields, if RF
%                   is not ON+OFF, only the subfield of the RF type is
%                   non-zero
%   .timeWeights    [ROIs x t_rf], optimal weights for fitted Gaussian
%                   spatial RF to match spatio-temporal RF (map)
%   .EV             [ROIs], explained variance of trace prediction based on
%                   Gaussian map (gaussMasks) * temporal weights
%                   (timeWeights)
%   .predictions    [t x ROIs], trace predictions based on Gaussian map
%                   (gaussMasks) * temporal weights (timeWeights)
%   .time_prediction [t], sample times of trace predictions
%   .x              [1 x rows], list of unique azimuth positions of circles
%   .y              [1 x columns], list of unique elevation positions of circles
%   .diameters      [1 x sizes], list of unique diameters of circles
%   .sizeTuning     [ROIs x sizes], size tuning curve
%   .time_RF        [t_rf], sample times of spatio-temporal RFs

data.maps = readNPY(fullfile(folder, 'circlesRf.maps.npy'));
data.bestSubFields = readNPY(fullfile(folder, 'circlesRf.bestSubField.npy'));
data.subFieldSigns = readNPY(fullfile(folder, 'circlesRf.subfieldSigns.npy'));
data.gaussPars = readNPY(fullfile(folder, 'circlesRf.gaussFitPars.npy'));
data.peaks = readNPY(fullfile(folder, 'circlesRf.peak.npy'));
data.gaussMasks = readNPY(fullfile(folder, 'circlesRf.gaussMask.npy'));
data.timeWeights = readNPY(fullfile(folder, 'circlesRf.gaussTimeWeights.npy'));
data.EV = readNPY(fullfile(folder, 'circlesRf.explVars.npy'));
data.predictions = readNPY(fullfile(folder, 'circlesRfPrediction.traces.npy'));
data.time_prediction = readNPY(fullfile(folder, 'circlesRfPrediction.timestamps.npy'));
data.x = readNPY(fullfile(folder, 'circlesRfDescr.x.npy'));
data.y = readNPY(fullfile(folder, 'circlesRfDescr.y.npy'));
data.diameters = readNPY(fullfile(folder, 'circlesRfDescr.diameters.npy'));
data.sizeTuning = readNPY(fullfile(folder, 'circlesRf.sizeTuning.npy'));
data.time_RF = readNPY(fullfile(folder, 'circlesRfDescr.timestamps.npy'));