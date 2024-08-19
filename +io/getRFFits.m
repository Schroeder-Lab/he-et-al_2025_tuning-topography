function data = getRFFits(f)
%GETRFFITS   Load receptive fields fitted to responses to visual noise.

% INPUTS
% folder            path to data of recording session

% OUTPUTS
% data
%   .maps           [ROIs x rows x columns x t_rf x type], spatio-temporal
%                   RF for ON (type=1) and OFF (type=2) subfields
%   .types          {ROIs}, 'ON', 'OFF', or 'ON+OFF'
%   .fitParameters  [ROIs x parameters], (amplitude, xCenter, xSTD,
%                   yCenter, ySTD, rotation)
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

data.maps = readNPY(fullfile(f, '_ss_rf.maps.npy'));
data.types = readcell(fullfile(f, '_ss_rf.type.csv'));
data.fitParameters = readNPY(fullfile(f, '_ss_rf.gaussFitPars.npy'));
data.peaks = readNPY(fullfile(f, '_ss_rf.peak.npy'));
data.gaussMasks = readNPY(fullfile(f, '_ss_rf.gaussMask.npy'));
data.timeWeights = readNPY(fullfile(f, '_ss_rf.gaussTimeWeights.npy'));
data.EV = readNPY(fullfile(f, '_ss_rf.explVars.npy'));
data.predictions = readNPY(fullfile(f, '_ss_rfPrediction.traces.npy'));
data.time_prediction = readNPY(fullfile(f, '_ss_rfPrediction.timestamps.npy'));
data.edges = readNPY(fullfile(f, '_ss_rfDescr.edges.npy'));
data.time_RF = readNPY(fullfile(f, '_ss_rfDescr.timestamps.npy'));