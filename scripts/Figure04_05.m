function Figure04_05(folders)

%% Parameters
sets = {'boutons', 'neurons'};
measures = {'direction', 'orientation'};
retinotopyRF = [false true]; % true: use RF positions estimated from 
                             % retinotopic mapping;
                             % false: use RF positions from RF mapping

%% For all plots
fPlots = fullfile(folders.plots, 'Figures', 'Figure04');
if ~isfolder(fPlots)
    mkdir(fPlots)
end

%% Load data: RF position, tuning preferences
% data: .rfPos, .oriPref, .OSI, .dirPref, .DSI, .set
data = Figure04_loadData(folders, sets, retinotopyRF);

%% Scatterplot showing preferred direction/orientation of each unit
Figure04_preferenceMapsAcrossAllDatasets(data, fPlots, sets, retinotopyRF)

%% Polar histograms of preferences per visual patch
Figure04_polarHists(data, fPlots, sets)

%% Plot smoothed preference maps pooling datasets
maps = Figure04_smoothedPreferenceMaps(data, fPlots, sets, retinotopyRF);

%% For all plots
fPlots = fullfile(folders.plots, 'Figures', 'Figure05');
if ~isfolder(fPlots)
    mkdir(fPlots)
end

%% Plot histograms + scatters: consistencies compared to null distribution
Figure05_consistenciesAll(maps, fPlots, sets, retinotopyRF, measures)

%% Plot consistencies compared to null distribution, per dataset
Figure05_consistenciesPerDataset(data, fPlots, sets, retinotopyRF, ...
    measures)