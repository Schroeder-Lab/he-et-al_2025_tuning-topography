function data = getRFFits(f)

data.maps = readNPY(fullfile(f, '_ss_rf.maps.npy'));
data.types = readcell(subFields, fullfile(f, '_ss_rf.type.csv'));
data.fitParameters = readNPY(fitPars, fullfile(f, '_ss_rf.gaussFitPars.npy'));
data.peaks = readNPY(peak_from_noise, fullfile(f, '_ss_rf.peak.npy'));
data.gaussMasks = readNPY(fitGaussians, fullfile(f, '_ss_rf.gaussMask.npy'));
data.timeWeights = readNPY(fitWeights, fullfile(f, '_ss_rf.gaussTimeWeights.npy'));
data.EV = readNPY(EVs, fullfile(f, '_ss_rf.explVars.npy'));
data.predictions = readNPY(predictions, fullfile(f, '_ss_rfPrediction.traces.npy'));
data.time_prediction = readNPY(t_pred, fullfile(f, '_ss_rfPrediction.timestamps.npy'));
data.edges = readNPY(stimPos, fullfile(f, '_ss_rfDescr.edges.npy'));
data.time_RF = readNPY(RFtimesInFrames * stimFrameDur, fullfile(f, '_ss_rfDescr.timestamps.npy'));