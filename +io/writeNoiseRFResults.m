function writeNoiseRFResults(results, folder)

writeNPY(results.maps, fullfile(folder, '_ss_rf.maps.npy'));
writeNPY(results.bestSubfields, fullfile(folder, '_ss_rf.bestSubField.npy'));
writeNPY(results.subfieldSigns, fullfile(folder, '_ss_rf.subfieldSigns.npy'));
writeNPY(results.gaussPars, fullfile(folder, '_ss_rf.gaussFitPars.npy'));
writeNPY(results.peaks, fullfile(folder, '_ss_rf.peak.npy'));
writeNPY(results.gaussMasks, fullfile(folder, '_ss_rf.gaussMask.npy'))
writeNPY(results.timeWeights, fullfile(folder, '_ss_rf.gaussTimeWeights.npy'))
writeNPY(results.EV, fullfile(folder, '_ss_rf.explVars.npy'));
writeNPY(results.units, fullfile(folder, '_ss_rf.clusters.npy'))
writeNPY(results.predictions, fullfile(folder, '_ss_rfPrediction.traces.npy'))
writeNPY(results.time_prediction, fullfile(folder, '_ss_rfPrediction.timestamps.npy'))
writeNPY(results.edges, fullfile(folder, '_ss_rfDescr.edges.npy'));
writeNPY(results.time_RF, fullfile(folder, '_ss_rfDescr.timestamps.npy'));