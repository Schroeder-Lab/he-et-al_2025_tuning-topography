function Figure03S(folders, glob)

%% Parameters
sets = {'boutons', 'neurons'};
retinotopyRF = [false true]; % true: use RF positions estimated from 
                             % retinotopic mapping;
                             % false: use RF positions from RF mapping
selectivityThresholds = [0.2 0.2; 0.1 0.2];

%% For all plots
fPlots = fullfile(folders.plots, 'Figures', 'Figure03S');
if ~isfolder(fPlots)
    mkdir(fPlots)
end

%% Load data: RF position, tuning preferences
% data: .rfPos, .oriPref, .OSI, .dirPref, .DSI, .set
data = Figures_loadData(folders, sets, retinotopyRF);

%% Polar histograms for DS and OS only units
Figure03_polarHists(glob, fPlots, data, sets, selectivityThresholds)