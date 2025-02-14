function Figure02(folders)

%% Parameters
sets = {'boutons', 'neurons'};
maxP = 0.05; % p-value threshold for response kernel and 
             % direction/orientation selectivity
minROIs = 15;
binSize = [5, 20];
stepSize = [2.5, 5];
xLims = [50 500];
cLims = [0.0006 0.00005];
fovLims = [20 160; 400 900];
numPerm = 1000;

%% Examples
ex = cell(2,2); % rows: (1) bouton, (2) neuron
ex(1,:) = {'SS078', '2017-10-05'};
% ex(1,:) = {'SS078', '2017-09-28', 1};
ex(2,:) = {'SS041', '2015-04-11'};
% ex(2,:) = {'SS044', '2015-04-28', 3};

%% For all plots
fPlots = fullfile(folders.plots, 'Figures', 'Figure02');
if ~isfolder(fPlots)
    mkdir(fPlots)
end

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

%% Plot pairwise distance in brain versus difference in tuning preference
for s = 1:2 % boutons and neurons
    subjDirs = dir(fullfile(folders.data, sets{s}, 'SS*'));
    dirDist = {};
    dirDiff = {};
    dirDiffNull = {};
    dirDiffRelative = {};
    distBinnedDir = {};
    oriDist = {};
    oriDiff = {};
    oriDiffNull = {};
    oriDiffRelative = {};
    distBinnedOri = {};
    fovSize = [];
    rec = 1;
    exRecs = [0 0];
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
            data = io.getRecordingInfo(f);
            roiPos = data.roiPositions(:,1:2);
            fovs = data.fovMicrons;
            [dirTuning, oriTuning] = io.getTuningResults(f, 'gratingsDrifting');

            dp = dirTuning.preference;
            validDir = ~isnan(dp) & dirTuning.pValue < maxP;
            if sum(validDir) < minROIs
                dirDist{rec} = [];
                dirDiff{rec} = [];
                dirDiffNull{rec} = [];
                dirDiffRelative{rec} = [];
                distBinnedDir{rec} = [];
            else
                % for all unit pairs, determine distance in brain (ignore
                % depth);
                ddist = spatial.determineDistance(roiPos(validDir,1), ...
                    roiPos(validDir,2));
                dp = dp(validDir);
                % for all unit pairs, determine difference between preferred
                % directions
                ddiff = tuning.determinePreferenceDiff(dp, 'dir');
                % permute preferences to test significance
                ddiffPermuted = NaN(length(ddiff), numPerm);
                rng('default');
                for k = 1:numPerm
                    order = randperm(length(dp));
                    ddiffPermuted(:,k) = tuning.determinePreferenceDiff(dp(order), 'dir');
                end
                % collect results
                dirDist{rec} = ddist;
                dirDiff{rec} = ddiff;
                dirDiffNull{rec} = ddiffPermuted;

                % difference in preference relative to null distribution
                [dirDiffRelative{rec}, distBinnedDir{rec}] = ...
                    spatial.getPrefDiffsRelativeNull(ddist, ddiff, ddiffPermuted, ...
                    binSize(s), stepSize(s), 'zscore');
            end
            op = oriTuning.preference;
            validOri = ~isnan(op) & oriTuning.pValue < maxP;
            if sum(validOri) < minROIs
                oriDist{rec} = [];
                oriDiff{rec} = [];
                oriDiffNull{rec} = [];
                oriDiffRelative{rec} = [];
                distBinnedOri{rec} = [];
            else
                % for all unit pairs, determine distance in brain (ignore
                % depth);
                odist = spatial.determineDistance(roiPos(validOri,1), ...
                    roiPos(validOri,2));
                op = op(validOri);
                % for all unit pairs, determine difference between preferred
                % orientations
                odiff = tuning.determinePreferenceDiff(op, 'ori');
                % permute preferences to test significance
                odiffPermuted = NaN(length(odiff), numPerm);
                rng('default');
                for k = 1:numPerm
                    order = randperm(length(op));
                    odiffPermuted(:,k) = tuning.determinePreferenceDiff(op(order), 'ori');
                end
                % collect results
                oriDist{rec} = odist;
                oriDiff{rec} = odiff;
                oriDiffNull{rec} = odiffPermuted;

                % difference in preference relative to null distribution
                [oriDiffRelative{rec}, distBinnedOri{rec}] = ...
                    spatial.getPrefDiffsRelativeNull(odist, odiff, odiffPermuted, ...
                    binSize(s), stepSize(s), 'zscore');
            end
            fovSize(rec) = mean(sqrt(sum(fovs.^2,2)));

            if strcmp(name,ex{s,1}) && strcmp(date,ex{s,2})
                exRecs(s) = rec;
            end
            rec = rec + 1;
        end
    end

    % plot distance vs tuning difference across all datasets
    n = sum(~any(isnan([cat(1, dirDist{:}) cat(1, dirDiff{:})]), 2));
    fig = spatial.plotPrefDiffVsDist(cat(1, dirDist{:}), ...
        cat(1, dirDiff{:}), cat(1, dirDiffNull{:}), ...
        binSize(s), stepSize(s), false);
    xlim([0 xLims(s)])
    ylim([0 180])
    clim([0 cLims(s)/2])
    title(['\DeltaDirection pref. vs \Deltaposition (n = ' num2str(n) ')'])
    io.saveFigure(fig, fPlots, ...
        sprintf('distanceAll_%s_direction', sets{s}))
    n = sum(~any(isnan([cat(1, oriDist{:}) cat(1, oriDiff{:})]), 2));
    fig = spatial.plotPrefDiffVsDist(cat(1, oriDist{:}), ...
        cat(1, oriDiff{:}), cat(1, oriDiffNull{:}), ...
        binSize(s), stepSize(s), false);
    xlim([0 xLims(s)])
    ylim([0 90])
    clim([0 cLims(s)])
    title(['\DeltaOrientation pref. vs \Deltaposition (n = ' num2str(n) ')'])
    io.saveFigure(fig, fPlots, ...
        sprintf('distanceAll_%s_orientation', sets{s}))

    % plot direction preference difference relative to null distribution (per
    % dataset)
    mini = ceil(min(cat(1, distBinnedDir{:})) / 2) * 2;
    maxi = floor(max(cat(1, distBinnedDir{:})) / 2) * 2;
    x = mini:0.5:maxi;
    y = NaN(length(x), length(dirDiffRelative));
    for rec = 1:length(dirDiffRelative)
        if all(isnan(dirDiffRelative{rec}))
            continue
        end
        ind1 = find(x >= min(distBinnedDir{rec}), 1);
        ind2 = find(x <= max(distBinnedDir{rec}), 1, 'last');
        y(ind1:ind2,rec) = interp1(distBinnedDir{rec}, ...
            dirDiffRelative{rec}, x(ind1:ind2), "pchip");
    end
    figure
    hold on
    fill([0 maxi maxi 0], [-3 -3 3 3], 'k', 'FaceColor', 'k', ...
        'FaceAlpha', 0.2, 'EdgeColor', 'none')
    plot([0 maxi],[0 0], 'k')
    p = plot(x, y);
    p(exRecs(s)).LineWidth = 2;
    legend(p,'Location','bestoutside')
    set(gca, "Box", "off", "ColorOrder", turbo(size(y,2)), ...
        "YTick", -12:3:12)
    plot(x, smoothdata(median(y,2,"omitnan"), "movmean", 5) , 'k', "LineWidth", 2)
    xlim([0 xLims(s)])
    ylim([-12 12])
    xlabel('Distance (um)')
    ylabel('\DeltaPreference (relative to null distribution)')
    title('\DeltaDirection pref. vs \Deltaposition')
    io.saveFigure(gcf, fPlots, ...
        sprintf('distancePerDataset_%s_direction', sets{s}))
    % plot orientation preference difference relative to null distribution (per
    % dataset)
    mini = ceil(min(cat(1, distBinnedOri{:})) / 2) * 2;
    maxi = floor(max(cat(1, distBinnedOri{:})) / 2) * 2;
    x = mini:0.5:maxi;
    y = NaN(length(x), length(oriDiffRelative));
    for rec = 1:length(oriDiffRelative)
        if all(isnan(oriDiffRelative{rec}))
            continue
        end
        ind1 = find(x >= min(distBinnedOri{rec}), 1);
        ind2 = find(x <= max(distBinnedOri{rec}), 1, 'last');
        y(ind1:ind2,rec) = interp1(distBinnedOri{rec}, ...
            oriDiffRelative{rec}, x(ind1:ind2), "pchip");
    end
    figure
    hold on
    fill([0 maxi maxi 0], [-3 -3 3 3], 'k', 'FaceColor', 'k', ...
        'FaceAlpha', 0.2, 'EdgeColor', 'none')
    plot([0 maxi],[0 0], 'k')
    p = plot(x, y);
    p(exRecs(s)).LineWidth = 2;
    legend(p,'Location','bestoutside')
    set(gca, "Box", "off", "ColorOrder", turbo(size(y,2)), ...
        "YTick", -12:3:12)
    plot(x, smoothdata(median(y,2,"omitnan"), "movmean", 5), 'k', "LineWidth", 2)
    xlim([0 xLims(s)])
    ylim([-12 12])
    xlabel('Distance (um)')
    ylabel('\DeltaPreference (relative to null distribution)')
    title('\DeltaOrientation pref. vs \Deltaposition')
    io.saveFigure(gcf, fPlots, ...
        sprintf('distancePerDataset_%s_orientation', sets{s}))

    % plot mean preference difference for each dataset against size of 
    % imaged field-of-view
    figure
    c = zeros(length(fovSize),3);
    c(exRecs(s),:) = [1 0 0];
    scatter(fovSize, cellfun(@mean, dirDiff, ...
        repmat({"omitnan"},1,length(dirDiff))), 36, c, 'filled');
    xlim(fovLims(s,:))
    ylim([0 90])
    xlabel('FOV diagonal (um)')
    ylabel('Mean \Deltadirection')
    title(sprintf('%s', sets{s}))
    io.saveFigure(gcf, fPlots, ...
        sprintf('prefDiffPerDataset_%s_direction', sets{s}))
    figure
    scatter(fovSize, cellfun(@mean, oriDiff, ...
        repmat({"omitnan"},1,length(oriDiff))), 36, c, 'filled')
    xlim(fovLims(s,:))
    ylim([0 45])
    xlabel('FOV diagonal (um)')
    ylabel('Mean \Deltaorientation')
    title(sprintf('%s', sets{s}))
    io.saveFigure(gcf, fPlots, ...
        sprintf('prefDiffPerDataset_%s_orientation', sets{s}))
end