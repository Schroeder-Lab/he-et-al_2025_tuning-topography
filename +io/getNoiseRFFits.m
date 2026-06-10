function data = getNoiseRFFits(folder, withEye)
%GETNOISERFFITS   Load receptive fields fitted to responses to visual noise.

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
%   .edges          [left rigth top bottom] of RF map
%   .time_RF        [t_rf], sample times of spatio-temporal RFs

if nargin < 2 || ~withEye
    prefixRF = '_ss_rf';
    prefixPred = '_ss_rfPrediction';
else
    prefixRF = '_ss_rf_eye';
    prefixPred = '_ss_rfPrediction_eye';

    data.eyePos = readNPY(fullfile(folder, '_ss_rfPrediction_eye.eyePos.npy'));
    data.A = readNPY(fullfile(folder, 'eyeToGaze.transformMatrix.npy'));
    data.defaultPos = readNPY(fullfile(folder, 'eyeToGaze.defaultEyePos.npy'));
end

data.maps = readNPY(fullfile(folder, ...
    [prefixRF '.maps.npy']));
data.bestSubFields = readNPY(fullfile(folder, ...
    [prefixRF '.bestSubField.npy']));
data.subFieldSigns = readNPY(fullfile(folder, ...
    [prefixRF '.subFieldSigns.npy']));
data.gaussPars = readNPY(fullfile(folder, ...
    [prefixRF '.gaussFitPars.npy']));
data.peaks = readNPY(fullfile(folder, ...
    [prefixRF '.peak.npy']));
data.gaussMasks = readNPY(fullfile(folder, ...
    [prefixRF '.gaussMask.npy']));
data.timeWeights = readNPY(fullfile(folder, ...
    [prefixRF '.gaussTimeWeights.npy']));
data.EV = readNPY(fullfile(folder, ...
    [prefixRF '.explVars.npy']));

data.predictions = readNPY(fullfile(folder, ...
    [prefixPred '.traces.npy']));
data.time_prediction = readNPY(fullfile(folder, ...
    [prefixPred '.timestamps.npy']));

data.edges = readNPY(fullfile(folder, '_ss_rfDescr.edges.npy'));
data.time_RF = readNPY(fullfile(folder, '_ss_rfDescr.timestamps.npy'));