function Figure06(folders, glob)

%% Parameters
maxP = 0.05; % p-value threshold

% to evaluate RFs
minEV = 0.01; % minimum explained variance
minPeak = 5; % minimum peak of RF (compared to noise)

%% Examples
% RFs
exAnimal = 'FG010';
exDate = '2024-10-19';
% exUnits = [375 361 291 235];
% depths:   24  76 160 209
exUnits = [361 491 282 559];
% depths:   76 118 158 250
exColors = lines(length(exUnits));

%% For all plots
fPlots = fullfile(folders.plots, 'Figures', 'Figure06');
if ~isfolder(fPlots)
    mkdir(fPlots)
end

%% Example spike waveshapes, firing traces and tuning curves

%% Load data: tuning preferences, RF position
% include all units within SC that are responsive to gratings
data = Figure06_loadData(folders, maxP, minEV, minPeak);

%% Plot tuning preferences and selectivity against depth in SC (per recording)
Figure06_plotTuningData(data, fPlots, glob)

%% Pairwise differences in tuning preferences and SC depth
Figure06_pairwiseDifferences(fPlots, glob, data)

%% Plot example RFs
exSession = strcmp(exAnimal, {data.animal}) & ...
    strcmp(exDate, {data.date});
Figure06_exampleRFs(folders, fPlots, glob, data(exSession), ...
    minEV, minPeak, exUnits, exColors)

%% Compare preferred direction / orientation to expected value based on RF location
Figure06_polarHistograms(data, fPlots)