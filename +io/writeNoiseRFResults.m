function writeNoiseRFResults(results, folder)
%WRITENOISERFRESULTS   Save results of RF fitting procedure.

% INPUTS
% results
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
% folder            path to data of recording session

writeNPY(results.maps, fullfile(folder, '_ss_rf.maps.npy'));
writeNPY(results.bestSubFields, fullfile(folder, '_ss_rf.bestSubField.npy'));
writeNPY(results.subFieldSigns, fullfile(folder, '_ss_rf.subfieldSigns.npy'));
writeNPY(results.gaussPars, fullfile(folder, '_ss_rf.gaussFitPars.npy'));
writeNPY(results.peaks, fullfile(folder, '_ss_rf.peak.npy'));
writeNPY(results.gaussMasks, fullfile(folder, '_ss_rf.gaussMask.npy'))
writeNPY(results.timeWeights, fullfile(folder, '_ss_rf.gaussTimeWeights.npy'))
writeNPY(results.EV, fullfile(folder, '_ss_rf.explVars.npy'))
writeNPY(results.predictions, fullfile(folder, '_ss_rfPrediction.traces.npy'))
writeNPY(results.time_prediction, fullfile(folder, '_ss_rfPrediction.timestamps.npy'))
writeNPY(results.edges, fullfile(folder, '_ss_rfDescr.edges.npy'));
writeNPY(results.time_RF, fullfile(folder, '_ss_rfDescr.timestamps.npy'));