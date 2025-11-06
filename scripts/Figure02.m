function Figure02(folders, glob)

%% Parameters
sets = {'boutons', 'neurons'};
maxP = 0.05; % p-value threshold for response kernel and 
             % direction/orientation selectivity

%% Examples
% RFs
ex = cell(2,4); % rows: (1) bouton, (2) neuron
% boutons:
% good retinotopic maps: SS077 2017-10-03, SS078_2017-10-05
ex(1,:) = {'SS078', '2017-10-05', 1, [9 62 146 77]};
% ex(1,:) = {'SS078', '2017-10-05', 1, [77 9 157 76 62 46 182 58 146 223 160 214 72]};
ex(2,:) = {'SS044', '2015-04-28', 3, [389 306 343 227]}; % also 427 (OFF)
% ex(2,:) = {'SS044', '2015-04-28', 3, [408 439 406 389]};
% ex(2,:) = {'SS044', '2015-04-28', 3, [368 385 395 286 383 408 278 390 354 387 367 369 326 306 389 321 348 378 421 252 422 268 295 294 233 308 439 242 432 285 267   327   361   406 343   403   427   227   416   380   420   335   392   372   386   313   376   277   282   417]};

% ROIs in imaging planes
% ex(1,:) = {'SS078', '2017-10-05'};
% % ex(1,:) = {'SS078', '2017-09-28', 1};
% % ex(2,:) = {'SS041', '2015-04-11'};
% ex(2,:) = {'SS044', '2015-04-28'};

retinotopyRF = [false true]; % true: use RF positions estimated from 
                             % retinotopic mapping;
                             % false: use RF positions from RF mapping

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

%% Example maps showing preferences of ROIs
for s = 1:2 % boutons and neurons
    str = sets{s};
    f = fullfile(folders.data, str, ex{s,1}, ex{s,2});
    % load data
    [dirTuning, oriTuning] = io.getTuningResults(f, 'gratingsDrifting');
    data = io.getRecordingInfo(f);
    masks = data.roiMasks;
    fovPix = data.fovPix;
    fovM = data.fovMicrons;

    tuning.plotOrientationMap(dirTuning.preference, ...
        dirTuning.pValue < maxP, 'dir', masks, fovPix(1,:), fovM(1,:));
    io.saveFigure(gcf, fPlots, sprintf('example_%s_directionMap_%s_%s', ...
        str, ex{s,1}, ex{s,2}))
    tuning.plotOrientationMap(oriTuning.preference, ...
        oriTuning.pValue < maxP, 'ori', masks, fovPix(1,:), fovM(1,:));
    io.saveFigure(gcf, fPlots, sprintf('example_%s_orientationMap_%s_%s', ...
        str, ex{s,1}, ex{s,2}))
end

%% Load data: RF position, tuning preferences
% data: .rfPos, .oriPref, .OSI, .dirPref, .DSI, .set
data = Figures_loadData(folders, sets, retinotopyRF);

%% Scatterplot showing preferred direction/orientation of each unit
Figure02_preferenceMapsAcrossAllDatasets(glob, fPlots, data, sets, ...
    retinotopyRF)