function Figure07(folders, glob)

%% Parameters
maxP = 0.05; % p-value threshold

% to evaluate RFs
minEV = 0.01; % minimum explained variance
minPeak = 5; % minimum peak of RF (compared to noise)

%% For all plots
fPlots = fullfile(folders.plots, 'Figures', 'Figure07');
if ~isfolder(fPlots)
    mkdir(fPlots)
end

%% Load data: tuning preferences, RF position
% include all units within SC that are responsive to gratings
data = Figure06_loadData(folders, maxP, minEV, minPeak);

%% Compare preferred direction / orientation to expected value based on RF location
Figure07_polarHistograms(glob, fPlots, data)