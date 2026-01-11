function writeCircleRFResults(results, folder)
%WRITECIRCLERFRESULTS   Save results of RF fitting procedure.

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
%   .x              [1 x rows], list of unique azimuth positions of circles
%   .y              [1 x columns], list of unique elevation positions of circles
%   .diameters      [1 x sizes], list of unique diameters of circles
%   .sizeTuning     [ROIs x sizes], size tuning curve
%   .time_RF        [t_rf], sample times of spatio-temporal RFs
% folder            path to data of recording session

writeNPY(results.maps, fullfile(folder, 'circlesRf.maps.npy'));
writeNPY(results.bestSubFields, fullfile(folder, 'circlesRf.bestSubField.npy'));
writeNPY(results.subFieldSigns, fullfile(folder, 'circlesRf.subfieldSigns.npy'));
writeNPY(results.gaussPars, fullfile(folder, 'circlesRf.gaussFitPars.npy'));
writeNPY(results.peaks, fullfile(folder, 'circlesRf.peak.npy'));
writeNPY(results.gaussMasks, fullfile(folder, 'circlesRf.gaussMask.npy'))
writeNPY(results.timeWeights, fullfile(folder, 'circlesRf.gaussTimeWeights.npy'))
writeNPY(results.EV, fullfile(folder, 'circlesRf.explVars.npy'))
writeNPY(results.predictions, fullfile(folder, 'circlesRfPrediction.traces.npy'))
writeNPY(results.time_prediction, fullfile(folder, 'circlesRfPrediction.timestamps.npy'))
writeNPY(results.x, fullfile(folder, 'circlesRfDescr.x.npy'));
writeNPY(results.y, fullfile(folder, 'circlesRfDescr.y.npy'));
writeNPY(results.diameters, fullfile(folder, 'circlesRfDescr.diameters.npy'));
writeNPY(results.sizeTuning, fullfile(folder, 'circlesRf.sizeTuning.npy'));
writeNPY(results.time_RF, fullfile(folder, 'circlesRfDescr.timestamps.npy'));