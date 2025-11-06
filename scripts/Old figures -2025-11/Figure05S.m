function Figure05S(folders)

%% Parameters
sets = {'boutons', 'neurons'};
retinotopyRF = [false true]; % true: use RF positions estimated from 
                             % retinotopic mapping;
                             % false: use RF positions from RF mapping
selectivityThresholds = [0.2 0.2; 0.1 0.2];

%% For all plots
fPlots = fullfile(folders.plots, 'Figures', 'Figure05S');
if ~isfolder(fPlots)
    mkdir(fPlots)
end

%% Load data: RF position, tuning preferences
% data: .rfPos, .oriPref, .OSI, .dirPref, .DSI, .set
data = Figure04_loadData(folders, sets, retinotopyRF);

%% Polar histograms for DS and OS only units
Figure05_polarHists(data, fPlots, sets, selectivityThresholds)