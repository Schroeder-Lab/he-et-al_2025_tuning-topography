function Figure03(folders)

%% Parameters
sets = {'boutons', 'neurons'};
maxP = 0.05; % p-value threshold for response kernel and 
             % direction/orientation selectivity
minROIs = 15;
binSize = [5, 20];
stepSize = [2.5, 5];
xLims = [50 500];
fovLims = [20 160; 400 900];
numPerm = 1000;

%% Examples
ex = cell(2,3); % rows: (1) bouton, (2) neuron
% boutons:
% good retinotopic maps: SS077 2017-10-03, SS078_2017-10-05
ex(1,:) = {'SS078', '2017-09-28', 1};
ex(2,:) = {'SS044', '2015-04-28', 3};

%% For all plots
fPlot = fullfile(folders.plots, 'Figure03');
if ~isfolder(fPlot)
    mkdir(fPlot)
end