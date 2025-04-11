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

%% Scatterplots for DS and OS only units
Figure04_preferenceMapsAcrossAllDatasets(data, fPlots, sets, ...
    retinotopyRF, selectivityThresholds)

%% Polar historgrams for DS and OS only units
Figure04_polarHists(data, fPlots, sets, selectivityThresholds)


%% Load retinotopic data for boutons
dataB = Figure04_loadData(folders, sets(1), true);

%% Polar histograms for direction and orientation preferences
Figure04_polarHists(dataB, fPlots, sets(1), [0 1; 0 1], '_retinotopic')