function Figure07(folders, glob)

%% Parameters
maxP = 0.05; % p-value threshold

% to evaluate RFs
minEV = 0.01; % minimum explained variance
minPeak = 5; % minimum peak of RF (compared to noise)

%% Examples
% RFs
exAnimal = 'FG010';
exDate = '2024-10-19';
exUnits = [361 491 282 559];
exColors = lines(length(exUnits));

%% For all plots
fPlots = fullfile(folders.plots, 'Figures', 'Figure07');
if ~isfolder(fPlots)
    mkdir(fPlots)
end

%% Load data: tuning preferences, RF position
% include all units within SC that are responsive to gratings
data = Figure06_loadData(folders, maxP, minEV, minPeak);

%% Plot example RFs
exSession = strcmp(exAnimal, {data.animal}) & ...
    strcmp(exDate, {data.date});
Figure07_exampleRFs(folders, fPlots, glob, data(exSession), ...
    minEV, minPeak, exUnits, exColors)

%% Compare preferred direction / orientation to expected value based on RF location
Figure07_polarHistograms(glob, fPlots, data)