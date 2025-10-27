function Figure06_pairwiseDifferences(fPlots, glob, data)

numPerm = 1000;
binSize = 50;
stepSize = 5;
xLims = 500;
cLims = 0.00003;
yTicks = {0:45:180, 0:45:90};
barWidth = 0.3;

features = {'direction', 'orientation'};
featNames = {'dir', 'ori'};

dist = cell(length(data), 2);
diff = cell(length(data), 2);
diffNull = cell(length(data), 2);
for feat = 1:2 % direction and orientation preferences
    for session = 1:length(data)
        units = data(session).([featNames{feat} 'Tuned']);
        % for all unit pairs, determine distance in depth;
        sdist = spatial.determineDistance( ...
            data(session).depth(units), zeros(sum(units),1));
        % for all unit pairs, determine difference between preferred
        % directions
        spref = data(session). ...
            ([featNames{feat} 'Preferences'])(units);
        sdiff = tuning.determinePreferenceDiff(spref, featNames{feat});
        % permute preferences to test significance
        sdiffPermuted = NaN(length(sdiff), numPerm);
        rng('default');
        for k = 1:numPerm
            order = randperm(sum(units));
            sdiffPermuted(:,k) = tuning.determinePreferenceDiff( ...
                spref(order), featNames{feat});
        end
        % collect results
        dist{session, feat} = sdist;
        diff{session, feat} = sdiff;
        diffNull{session, feat} = sdiffPermuted;
    end
end

%% Plot distance vs tuning difference across all datasets
for feat = 1:2
    n = sum(~any(isnan([cat(1, dist{:, feat}) cat(1, diff{:, feat})]), 2));
    fig = spatial.plotPrefDiffVsDist(cat(1, dist{:, feat}), ...
        cat(1, diff{:, feat}), cat(1, diffNull{:, feat}), ...
        binSize, stepSize, false);
    set(gcf, 'Position', glob.figPositionDefault)
    set(gca, 'YTick', yTicks{feat})
    xlim([0 xLims])
    ylim([0 180/feat])
    clim([0 cLims/(3-feat)])
    title(['\Delta' features{feat} ' pref. vs \Deltaposition (n = ' ...
        num2str(n) ')'])
    io.saveFigure(fig, fPlots, sprintf('distanceAll_%s', features{feat}))
end

%% Plot mean pairwise preference difference in each recording
% (and permutation test with preference permutation across recordings)
[~, ~, mouseID] = unique({data.animal});
tuned = [cat(1, data.dirTuned), cat(1, data.oriTuned)];
preferences = [cat(1, data.dirPreferences), ...
    cat(1, data.oriPreferences)];

for feat = 1:2
    meanDiffs = cellfun(@mean, diff(:,feat));
    meanDiffs_shuffled = NaN(length(data), numPerm);
    validUnits = find(tuned(:, feat));
    prefs = preferences(validUnits, feat);
    for session = 1:length(data)
        rng('default');
        for k = 1:numPerm
            indPerm = randsample(length(validUnits), ...
                sum(data(session).([featNames{feat} 'Tuned'])));
            diffPerm = tuning.determinePreferenceDiff(prefs(indPerm), ...
                featNames{feat});
            meanDiffs_shuffled(session, k) = mean(diffPerm);
        end
    end
    confIntv = prctile(meanDiffs_shuffled, [2.5 97.5], 2);
    medians = median(meanDiffs_shuffled, 2);
    signif = find(meanDiffs < confIntv(:,1) | meanDiffs > confIntv(:,2));
    figure('Position', [100 100 1090 300])
    hold on
    for session = 1:length(data)
        fill(session + [-1 1 1 -1].* barWidth, ...
            reshape(repmat(confIntv(session,:), 2, 1), [], 1), 'k', ...
            'EdgeColor', 'none', 'FaceColor', [1 1 1] .* 0.8)
        plot(session + [-1 1] .* barWidth, [1 1] .* medians(session), ...
            'Color', [1 1 1] .* 0.5, 'LineWidth',1)
    end
    scatter(1:length(data), meanDiffs, 'k', 'filled')
    scatter(signif, ones(size(signif)) .* (180 / feat), 30, 'k', '*')
    xlim([0 length(data) + 1])
    ylim([0 180 / feat])
    set(gca, "XTick", 1:length(data), "XTickLabel", mouseID, ...
        "YTick", 0 : (30 / feat) : (180 / feat))
    xlabel('Animal ID')
    ylabel(['\Delta' features{feat} ' preference'])
    io.saveFigure(gcf, fPlots, ...
        sprintf('tuningDifference_%s', features{feat}))
end