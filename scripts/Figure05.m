function Figure05(folders, glob)

%% Parameters
sets = {'boutons', 'neurons'};
% for evaluation of receptive fields
retinotopyRF = [false true]; % true: use RF positions estimated from 
                             % retinotopic mapping;
                             % false: use RF positions from RF mapping
minEV = 0.01; % minimum explained variance to plot RF
minPeak = 5; % minimum peak of RF (compared to noise) to plot RF
dist2edge = 5; % minimum distance of RF centre to monitor edge

fovLims = [20 160; 400 900];
maxP = 0.05; % p-value threshold for direction/orientation selectivity
minROIs = 15;
numPerm = 1000;
binSize = [5, 20];
stepSize = [2.5, 5];
xLims = [50 400];
cLims = [0.0004 0.000044];

%% Examples
% ROIs in imaging planes
ex(1,:) = {'SS078', '2017-10-05'};
ex(2,:) = {'SS041', '2015-04-11'};

%% For all plots
fPlots = fullfile(folders.plots, 'Figures', 'Figure05');
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
            if ~isfile(fullfile(f, '_ss_gratingsDrifting.intervals.npy')) || ...
                    ~isfile(fullfile(f, '_ss_sparseNoise.times.npy'))
                continue
            end
                
            % load data
            data = io.getRecordingInfo(f);
            roiPos = data.roiPositions(:,1:2);
            fovs = data.fovMicrons;
            [dirTuning, oriTuning] = io.getTuningResults(f, 'gratingsDrifting');
            if ~retinotopyRF(s)
                dataRF = io.getNoiseRFFits(f);
                edges = dataRF.edges; % [left right top bottom]
                pos = dataRF.gaussPars(:,[2 4]);
                % exclude all RFs too close to monitor edge
                invalidRF = dataRF.EV < minEV | dataRF.peaks < minPeak | ...
                    pos(:,1) < edges(1)+dist2edge | ...
                    pos(:,1) > edges(2)-dist2edge | ...
                    pos(:,2) > edges(3)-dist2edge | ...
                    pos(:,2) < edges(4)+dist2edge;
                clear dataRF
            else
                if ~isfile(fullfile(f, '_ss_rf.posRetinotopy.npy'))
                    continue
                end
                pos = readNPY(fullfile(f, '_ss_rf.posRetinotopy.npy'));
                edges = readNPY(fullfile(f, '_ss_rfDescr.edges.npy'));
                % exclude all RFs too close to monitor edge
                invalidRF = pos(:,1) < edges(1)+dist2edge | ...
                    pos(:,1) > edges(2)-dist2edge | ...
                    pos(:,2) > edges(3)-dist2edge | ...
                    pos(:,2) < edges(4)+dist2edge;
            end

            dp = dirTuning.preference;
            validDir = ~isnan(dp) & dirTuning.pValue < maxP;
            valid = validDir & ~invalidRF;
            if sum(valid) < minROIs
                dirDist{rec} = [];
                dirDiff{rec} = [];
                dirDiffNull{rec} = [];
                dirDiffRelative{rec} = [];
                distBinnedDir{rec} = [];
            else
                % for all unit pairs, determine distance in brain (ignore
                % depth);
                ddist = spatial.determineDistance(roiPos(valid,1), ...
                    roiPos(valid,2));
                dp = dp(valid);
                % for all unit pairs, determine difference between preferred
                % directions
                ddiff = tuning.determinePreferenceDiff(dp, 'dir');
                % permute preferences to test significance
                ddiffPermuted = NaN(length(ddiff), numPerm);
                rng('default');
                for k = 1:numPerm
                    order = randperm(length(dp));
                    ddiffPermuted(:,k) = ...
                        tuning.determinePreferenceDiff(dp(order), 'dir');
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
            valid = validOri & ~invalidRF;
            if sum(valid) < minROIs
                oriDist{rec} = [];
                oriDiff{rec} = [];
                oriDiffNull{rec} = [];
                oriDiffRelative{rec} = [];
                distBinnedOri{rec} = [];
            else
                % for all unit pairs, determine distance in brain (ignore
                % depth);
                odist = spatial.determineDistance(roiPos(valid,1), ...
                    roiPos(valid,2));
                op = op(valid);
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
    fig.Position = glob.figPositionDefault;
    set(gca, 'YTick', 0:45:180)
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
    fig.Position = glob.figPositionDefault;
    set(gca, 'YTick', 0:45:90)
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
    
    figure('Position', glob.figPositionDefault)
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
    title(['\DeltaDirection pref. vs \Deltaposition (n = ' ...
        num2str(sum(~all(isnan(y),1))) ')'])
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

    figure('Position', glob.figPositionDefault)
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
    title(['\DeltaOrientation pref. vs \Deltaposition (n = ' ...
        num2str(sum(~all(isnan(y),1))) ')'])
    io.saveFigure(gcf, fPlots, ...
        sprintf('distancePerDataset_%s_orientation', sets{s}))

    % plot mean preference difference for each dataset against size of 
    % imaged field-of-view
    figure('Position', glob.figPositionDefault)
    c = zeros(length(fovSize),3);
    c(exRecs(s),:) = [1 0 0];
    m = cellfun(@mean, dirDiff, repmat({"omitnan"},1,length(dirDiff)));
    scatter(fovSize, m, 36, c, 'filled');
    xlim(fovLims(s,:))
    ylim([0 100])
    xlabel('FOV diagonal (um)')
    ylabel('Mean \Deltadirection')
    title(sprintf('%s (n = %d)', sets{s}, sum(~isnan(m))))
    io.saveFigure(gcf, fPlots, ...
        sprintf('prefDiffPerDataset_%s_direction', sets{s}))

    figure('Position', glob.figPositionDefault)
    m = cellfun(@mean, oriDiff, repmat({"omitnan"},1,length(oriDiff)));
    scatter(fovSize, m, 36, c, 'filled')
    xlim(fovLims(s,:))
    ylim([0 45])
    xlabel('FOV diagonal (um)')
    ylabel('Mean \Deltaorientation')
    title(sprintf('%s (n = %d)', sets{s}, sum(~isnan(m))))
    io.saveFigure(gcf, fPlots, ...
        sprintf('prefDiffPerDataset_%s_orientation', sets{s}))
end