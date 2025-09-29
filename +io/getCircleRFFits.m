function data = getCircleRFFits(folder)

data.maps = readNPY(fullfile(folder, 'circlesRf.maps.npy'));
data.bestSubFields = readNPY(fullfile(folder, 'circlesRf.bestSubField.npy'));
data.subfieldSigns = readNPY(fullfile(folder, 'circlesRf.subfieldSigns.npy'));
data.fitParameters = readNPY(fullfile(folder, 'circlesRf.gaussFitPars.npy'));
data.peaks = readNPY(fullfile(folder, 'circlesRf.peak.npy'));
data.gaussMasks = readNPY(fullfile(folder, 'circlesRf.gaussMask.npy'));
data.timeWeights = readNPY(fullfile(folder, 'circlesRf.gaussTimeWeights.npy'));
data.EV = readNPY(fullfile(folder, 'circlesRf.explVars.npy'));
data.units = readNPY(fullfile(folder, 'circlesRf.clusters.npy'));
data.predictions = readNPY(fullfile(folder, 'circlesRfPrediction.traces.npy'));
data.time_prediction = readNPY(fullfile(folder, 'circlesRfPrediction.timestamps.npy'));
data.x = readNPY(fullfile(folder, 'circlesRfDescr.x.npy'));
data.y = readNPY(fullfile(folder, 'circlesRfDescr.y.npy'));
data.diameters = readNPY(fullfile(folder, 'circlesRfDescr.diameters.npy'));
data.sizeTuning = readNPY(fullfile(folder, 'circlesRf.sizeTuning.npy'));
data.time_RF = readNPY(fullfile(folder, 'circlesRfDescr.timestamps.npy'));