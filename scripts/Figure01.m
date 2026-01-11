function Figure01(folders, glob)

%% Parameters
sets = {'boutons', 'neurons'};
maxP = 0.05; % p-value threshold for response kernel and 
             % direction/orientation selectivity
minR2 = 0; % threshold for explained variance for response kernel

%% Examples
ex = cell(2,4); % rows: (1) bouton, (2) neuron
ex(1,:) = {'SS078', '2017-09-28', 1, [40 48  51 62 192  42]};
ex(2,:) = {'SS044', '2015-04-28', 3, [258 227 256 281 328 369]};

%% For all plots
fPlots = fullfile(folders.plots, 'Figures', 'Figure01');
if ~isfolder(fPlots)
    mkdir(fPlots)
end

%% Example mean frames, calcium traces and tuning curves
Figure01_examples(folders, fPlots, sets, ex);

%% Population direction tuning curves
Figure01_tuning_prefs_selectivity(folders, glob, fPlots, sets, maxP, ...
    minR2, ex);