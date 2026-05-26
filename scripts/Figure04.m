function Figure05(folders, glob)

%% Parameters
sets = {'boutons', 'neurons'};
measures = {'direction', 'orientation'};
retinotopyRF = [false true]; % true: use RF positions estimated from 
                             % retinotopic mapping;
                             % false: use RF positions from RF mapping

%% For all plots
fPlots = fullfile(folders.plots, 'Figures', 'Figure05');
if ~isfolder(fPlots)
    mkdir(fPlots)
end

%% Load data: RF position, tuning preferences
% data: .rfPos, .oriPref, .OSI, .dirPref, .DSI, .set
data = Figures_loadData(folders, sets, retinotopyRF);

%% Scatterplot showing preferred direction/orientation of each unit
Figure05_preferenceMapsAcrossAllDatasets(glob, fPlots, data, sets, ...
    retinotopyRF)

%% Plot smoothed preference maps pooling datasets
maps = Figure05_smoothedPreferenceMaps(glob, fPlots, data, sets, ...
    retinotopyRF);

%% Plot histograms: consistencies compared to null distribution
Figure05_consistenciesAll(glob, fPlots, maps, sets, retinotopyRF, measures)

%% Plot consistencies compared to null distribution, per dataset
Figure05_consistenciesPerDataset(glob, fPlots, data, sets, ...
    retinotopyRF, measures)

