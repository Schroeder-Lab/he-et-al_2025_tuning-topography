function Figure06(folders, glob)

%% Parameters
% for evaluation of tuning selectivity
maxP = 0.05;
minSI = 0.1;

% to evaluate RFs
minEV = 0.01; % minimum explained variance
minPeak = 5; % minimum peak of RF (compared to noise)
dist2edge = 5; % minimum distance of RF centre to monitor edge

%% Examples
% RFs
exAnimal = 'FG010';
exDate = '2024-10-19';
exUnits = [361 491 235];

%% For all plots
fPlots = fullfile(folders.plots, 'Figures', 'Figure06');
if ~isfolder(fPlots)
    mkdir(fPlots)
end

%% Example firing traces and tuning curves
Figure06_examplesTuning(folders, fPlots, exAnimal, exDate, exUnits)

%% Load data: tuning preferences, RF position
% include all units within SC that are responsive to gratings
data = Figure06_loadData(folders, maxP, minSI, minEV, minPeak, dist2edge);

%% Plot tuning preferences and selectivity against depth in SC (per recording)
Figure06_plotTuningData(fPlots, glob, data)

%% Pairwise differences in tuning preferences and SC depth
Figure06_pairwiseDifferences(fPlots, glob, data)