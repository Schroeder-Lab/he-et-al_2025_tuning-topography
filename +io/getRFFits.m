function data = getRFFits(f)

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