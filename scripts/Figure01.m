%% Folders
getFolders;

%% Parameters
sets = {'boutons', 'neurons'};
maxP = 0.05; % p-value threshold for response kernel and 
             % direction/orientation selectivity
minR2 = 0.02; % threshold for explained variance for response kernel

%% Examples
ex = cell(2,4); % rows: (1) bouton, (2) neuron
% final selection (6 neighbours)
ex(1,:) = {'SS078', '2017-09-28', 1, [40 48  51 62 192  42]};
% final selection (2 triples of neighbours)
ex(2,:) = {'SS044', '2015-04-28', 3, [258 227 256 281 328 369]};
% ex(2,:) = {'SS044', '2015-04-28', 3, [328 378 369 236 245 258]};
% good boutons
% ex(1,:) = {'SS078', '2017-09-28', 1, [1 3 4 5 6 8 9 10 12 14 15 16 17 18 19 20 21 23 24 26 27 29 31 32 33 34 35 36 38 39 40 42 45 46 48 49 50 51 52 57 58 60 62 63 64 71 72 74 80 83 95 99 100 105 107 109 110 123 125 126 128 135 142 145 146 150 165 168 175 176 180 181 186 187 192 196 200 201 205 207 208 217 224 229 239 242 252 257 262 269 276]};
% (almost)final selection
% ex(1,:) = {'SS078', '2017-09-28', 1, [40 46 48  51 62 192  42  74 205]};
% best OS boutons
% ex(1,:) = {'SS078', '2017-09-28', 1, [18 45 150]};
% other good datasets
% ex(1,:) = {'SS078', '2017-10-05', 1, []};
% ex(1,:) = {'SS077', '2017-10-03', 1, []};
% good neurons
% ex(2,:) = {'SS044', '2015-04-28', 3, [227 231 236 240 245 252 256 258 268 281 285 293 313 321 328 369 378 379 380 383 385 393 408 425 426 430]};
% best OS & DS neurons 231 245 293 321 378 
% ex(2,:) = {'SS044', '2015-04-28', 3, [227 240 252 256 258 268 281 313 328 369 380 385 393]};
% best OS neurons
% ex(2,:) = {'SS044', '2015-04-28', 3, [227 240]};

%% Add paths
addpath(genpath(fullfile(folders.tools, 'npy-matlab')))
addpath(fullfile(folders.repo))

%% For all plots
fPlot = fullfile(folders.plots, 'Figure01');
if ~isfolder(fPlot)
    mkdir(fPlot)
end

%% Example mean frames, calcium traces and tuning curves
Figure01_examples;

%% Population direction tuning curves
Figure01_tuning_prefs_selectivity;