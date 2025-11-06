function Figure04S(folders)

%% Parameters
sets = {'boutons', 'neurons'};
retinotopyRF = [false true]; % true: use RF positions estimated from 
                             % retinotopic mapping;
                             % false: use RF positions from RF mapping
selectivityThresholds = [0.2 0.2; 0.1 0.2];

%% For all plots
fPlots = fullfile(folders.plots, 'Figures', 'Figure04S');
if ~isfolder(fPlots)
    mkdir(fPlots)
end

%% Load data: RF position, tuning preferences
% data: .rfPos, .oriPref, .OSI, .dirPref, .DSI, .set
data = Figure04_loadData(folders, sets, retinotopyRF);

%% Scatterplots of preferences for DS and OS only units
Figure04_preferenceMapsAcrossAllDatasets(data, fPlots, sets, ...
    retinotopyRF, selectivityThresholds)

%% Plot smoothed preference maps pooling datasets
maps = Figure04_smoothedPreferenceMaps(data, [], sets, retinotopyRF);

%% Scatterplots: consistencies compared to null distribution
Figure04_consistenciesAll(maps, fPlots, sets, retinotopyRF, measures)