%% Folders
getFolders;

%% Parameters
sets = {'boutons', 'neurons'};
maxP = 0.05; % p-value threshold for response kernel and 
             % direction/orientation selectivity
minROIs = 15;
binSize = [10, 5];
stepSize = [5, 2.5];
numPerm = 1000;

%% Examples
ex = cell(2,3); % rows: (1) bouton, (2) neuron
ex(1,:) = {'SS078', '2017-09-28', 1};
ex(2,:) = {'SS044', '2015-04-28', 3};

%% Add paths
addpath(genpath(fullfile(folders.tools, 'npy-matlab')))
addpath(fullfile(folders.repo))

%% For all plots
fPlot = fullfile(folders.plots, 'Figure02');
if ~isfolder(fPlot)
    mkdir(fPlot)
end

%% Example maps showing preferences of ROIs
% for s = 1:2 % boutons and neurons
%     str = sets{s};
%     f = fullfile(folders.data, str, ex{s,1}, ex{s,2});
%     % load data
%     [dirTuning, oriTuning] = io.getTuningResults(f, 'gratingsDrifting');
%     data = io.getCalciumData(f);
%     planes = data.planes;
%     data = io.getRecordingInfo(f);
%     masks = data.roiMasks;
%     fovPix = data.fovPix;
%     fovM = data.fovMicrons;
% 
%     indP = planes == ex{s,3};
%     tuning.plotOrientationMap(dirTuning.preference(indP), ...
%         dirTuning.pValue(indP) < maxP, 'dir', masks(indP,:), ...
%         fovPix(ex{s,3},:), fovM(ex{s,3},:));
%     io.saveFigure(gcf, fPlot, sprintf('example_%s_directionMap_%s_%s_plane%02d', ...
%         str, ex{s,1}, ex{s,2}, ex{s,3}))
%     tuning.plotOrientationMap(oriTuning.preference(indP), ...
%         oriTuning.pValue(indP) < maxP, 'ori', masks(indP,:), ...
%         fovPix(ex{s,3},:), fovM(ex{s,3},:));
%     io.saveFigure(gcf, fPlot, sprintf('example_%s_orientationMap_%s_%s_plane%02d', ...
%         str, ex{s,1}, ex{s,2}, ex{s,3}))
% end

%% Plot pairwise distance in brain versus difference in tuning preference
for s = 1:2 % boutons and neurons
    subjDirs = dir(fullfile(folders.data, sets{s}, 'SS*'));
    dist = {};
    dirDiff = {};
    dirDiffNull = {};
    dirDiffRelative = {};
    distBinnedDir = {};
    oriDiff = {};
    oriDiffNull = {};
    oriDiffRelative = {};
    distBinnedOri = {};
    rec = 1;
    for subj = 1:length(subjDirs) % animals
        name = subjDirs(subj).name;
        fprintf('%s\n', name)
        dateDirs = dir(fullfile(folders.data, sets{s}, name, '2*'));
        for dt = 1:length(dateDirs) %dates
            date = dateDirs(dt).name;
            f = fullfile(folders.data, sets{s}, name, date);
            % ignore session if stimulus was not presented
            if ~isfile(fullfile(f, '_ss_gratingsDrifting.intervals.npy'))
                continue
            end
                
            % load data
            data = io.getCalciumData(f);
            planes = data.planes;
            data = io.getRecordingInfo(f);
            roiPos = data.roiPositions(:,1:2);
            [dirTuning, oriTuning] = io.getTuningResults(f, 'gratingsDrifting');

            dp = dirTuning.preference;
            op = oriTuning.preference;
            validDir = ~isnan(dp) & dirTuning.pValue < maxP;
            validOri = ~isnan(op) & oriTuning.pValue < maxP;
            indValid = validDir | validOri;
            if sum(indValid) < minROIs
                continue
            end
            if sum(validDir) < minROIs
                dp = NaN;
            else
                dp = dp(indValid);
            end
            if sum(validOri) < minROIs
                op = NaN;
            else
                op = op(indValid);
            end

            % for all unit pairs, determine difference between preferred
            % directions/orientations
            dd = tuning.determinePreferenceDiff(dp, 'dir');
            od = tuning.determinePreferenceDiff(op, 'ori');
            % for all unit pairs, determine distance in brain (ignore
            % depth);
            d = spatial.determineDistance(roiPos(indValid,1), ...
                roiPos(indValid,2));
            % permute preferences to test significance
            ddPermuted = NaN(length(dd), numPerm);
            odPermuted = NaN(length(od), numPerm);
            rng('default');
            for k = 1:numPerm
                order = randperm(length(dp));
                ddPermuted(:,k) = tuning.determinePreferenceDiff(dp(order), 'dir');
                odPermuted(:,k) = tuning.determinePreferenceDiff(op(order), 'ori');
            end
            % collect results
            dist{rec} = d;
            dirDiff{rec} = dd;
            oriDiff{rec} = od;
            dirDiffNull{rec} = ddPermuted;
            oriDiffNull{rec} = odPermuted;

            % difference in preference relative to null distribution
            [dirDiffRelative{rec}, distBinnedDir{rec}] = ...
                spatial.getPrefDiffsRelativeNull(d, dd, ddPermuted, ...
                binSize(s), stepSize(s));
            [oriDiffRelative{rec}, distBinnedOri{rec}] = ...
                spatial.getPrefDiffsRelativeNull(d, od, odPermuted, ...
                binSize(s), stepSize(s));

            rec = rec + 1;
        end
    end

    % plot across all datasets
    fig = spatial.plotPrefDiffVsDist(cat(1, dist{:}), ...
        cat(1, dirDiff{:}), cat(1, dirDiffNull{:}), ...
        binSize(s), stepSize(s), false);
    title('\DeltaDirection pref. vs \Deltaposition')
    io.saveFigure(fig, fPlot, ...
        sprintf('distanceAll_%s_direction', sets{s}))
    fig = spatial.plotPrefDiffVsDist(cat(1, dist{:}), ...
        cat(1, oriDiff{:}), cat(1, oriDiffNull{:}), ...
        binSize(s), stepSize(s), false);
    title('\DeltaOrientation pref. vs \Deltaposition')
    io.saveFigure(fig, fPlot, ...
        sprintf('distanceAll_%s_orientation', sets{s}))

    % plot direction preference difference relative to null distribution (per
    % dataset)
    mini = ceil(min(cat(1, distBinnedDir{:})) / 2) * 2;
    maxi = floor(max(cat(1, distBinnedDir{:})) / 2) * 2;
    x = mini:0.5:maxi;
    y = NaN(length(x), length(dirDiffRelative));
    for rec = 1:length(dirDiffRelative)
        ind1 = find(x >= min(distBinnedDir{rec}), 1);
        ind2 = find(x <= max(distBinnedDir{rec}), 1, 'last');
        y(ind1:ind2,rec) = interp1(distBinnedDir{rec}, ...
            dirDiffRelative{rec}, x(ind1:ind2), "pchip");
    end
    figure
    plot(x, y)
    set(gca, "Box", "off", "ColorOrder", turbo(size(y,2)))
    hold on
    plot(x, median(y,2,"omitnan"), 'k', "LineWidth", 2)
    plot([0 maxi],[1 1].*0.025, 'k')
    plot([0 maxi],[1 1].*0.975, 'k')
    xlim([0 maxi])
    xlabel('Distance (um)')
    ylabel('\DeltaPreference (relative to null distribution)')
    title('\DeltaDirection pref. vs \Deltaposition')
    io.saveFigure(gcf, fPlot, ...
        sprintf('distancePerDataset_%s_direction', sets{s}))
    % plot orientation preference difference relative to null distribution (per
    % dataset)
    mini = ceil(min(cat(1, distBinnedDir{:})) / 2) * 2;
    maxi = floor(max(cat(1, distBinnedDir{:})) / 2) * 2;
    x = mini:0.5:maxi;
    y = NaN(length(x), length(dirDiffRelative));
    for rec = 1:length(dirDiffRelative)
        ind1 = find(x >= min(distBinnedDir{rec}), 1);
        ind2 = find(x <= max(distBinnedDir{rec}), 1, 'last');
        y(ind1:ind2,rec) = interp1(distBinnedDir{rec}, ...
            dirDiffRelative{rec}, x(ind1:ind2), "pchip");
    end
    figure
    plot(x, y)
    set(gca, "Box", "off", "ColorOrder", turbo(size(y,2)))
    hold on
    plot(x, median(y,2,"omitnan"), 'k', "LineWidth", 2)
    plot([0 maxi],[1 1].*0.025, 'k')
    plot([0 maxi],[1 1].*0.975, 'k')
    xlim([0 maxi])
    xlabel('Distance (um)')
    ylabel('\DeltaPreference (relative to null distribution)')
    title('\DeltaOrientation pref. vs \Deltaposition')
    io.saveFigure(gcf, fPlot, ...
        sprintf('distancePerDataset_%s_orientation', sets{s}))
end