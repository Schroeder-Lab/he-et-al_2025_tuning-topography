function Figure02(folders, glob)

%% Parameters
sets = {'boutons', 'neurons'};

%% Examples
% RFs
ex = cell(2,4); % rows: (1) bouton, (2) neuron
ex(1,:) = {'SS078', '2017-10-05', 1, [9 62 146]};
ex(2,:) = {'SS044', '2015-04-28', 3, [389 306 227]};

%% For all plots
fPlots = fullfile(folders.plots, 'Figures', 'Figure02');
if ~isfolder(fPlots)
    mkdir(fPlots)
end

%% Examples: mean 2P images, RFs, RF outlines
Figure02_RFs(folders, sets, ex, fPlots)

%% Histograms: 
% Distance between measured and fitted retinotopic RF positions +
% Sizes of fitted RFs
Figure02_measured_vs_retinotopic_RF(folders, glob, sets, fPlots)