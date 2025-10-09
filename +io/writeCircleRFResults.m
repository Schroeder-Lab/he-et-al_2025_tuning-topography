function writeCircleRFResults(results, folder)

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