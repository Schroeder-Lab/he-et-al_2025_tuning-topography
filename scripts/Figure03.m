function Figure03(folders)

%% Parameters
sets = {'boutons', 'neurons'};

%% Examples
ex = cell(2,4); % rows: (1) bouton, (2) neuron
% boutons:
% good retinotopic maps: SS077 2017-10-03, SS078_2017-10-05
ex(1,:) = {'SS078', '2017-10-05', 1, [9 62 146 77]};
% ex(1,:) = {'SS078', '2017-10-05', 1, [77 9 157 76 62 46 182 58 146 223 160 214 72]};
ex(2,:) = {'SS044', '2015-04-28', 3, [389 306 343 227]}; % also 427 (OFF)
% ex(2,:) = {'SS044', '2015-04-28', 3, [408 439 406 389]};
% ex(2,:) = {'SS044', '2015-04-28', 3, [368 385 395 286 383 408 278 390 354 387 367 369 326 306 389 321 348 378 421 252 422 268 295 294 233 308 439 242 432 285 267   327   361   406 343   403   427   227   416   380   420   335   392   372   386   313   376   277   282   417]};

%% For all plots
fPlots = fullfile(folders.plots, 'Figures', 'Figure03');
if ~isfolder(fPlots)
    mkdir(fPlots)
end

%% Examples: mean 2P images, RFs, RF outlines
Figure03_examples(folders, sets, ex, fPlots)

%% Histograms: 
% Distance between measured and fitted retinotopic RF positions +
% Sizes of fitted RFs
Figure03_measured_vs_retinotopic_RF(folders, sets, fPlots)

%% Scatterplots: 
% Pairwise difference in pref. dir./ori. dependent on distance between RFs
Figure03_prefDiff_vs_distance(folders, sets, fPlots)