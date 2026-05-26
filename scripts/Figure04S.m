function Figure05S(folders, glob)

%% Parameters
sets = {'boutons', 'neurons'};
measures = {'direction', 'orientation'};
retinotopyRF = [false true]; % true: use RF positions estimated from 
                             % retinotopic mapping;
                             % false: use RF positions from RF mapping

%% For all plots
fPlots = fullfile(folders.plots, 'Figures', 'Figure05S');
if ~isfolder(fPlots)
    mkdir(fPlots)
end

%% Load data: RF position, tuning preferences
% data: .rfPos, .oriPref, .OSI, .dirPref, .DSI, .set
data = Figures_loadData(folders, sets, retinotopyRF);

%% Plot smoothed preference maps pooling datasets
maps = Figure05_smoothedPreferenceMaps(glob, [], data, sets, retinotopyRF);

%% Scatterplots: consistencies compared to null distribution
Figure05S_consistenciesAll(glob, fPlots, maps, sets, retinotopyRF, ...
    measures)