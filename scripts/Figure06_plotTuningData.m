function Figure06_plotTuningData(tuningData, fPlots)

%% Parameters
% tuning curves
dirBins = 0:30:360;
dirEdges = (0:30:390) - 15;
dirBinsFine = 0:1:360;
oriBins = 0:15:180;
oriEdges = (0:15:195) - 7.5;
oriBinsFine = 0:1:180;
minUnits = 10;

% preferences across depths
scale = 30;
colors = {colmaps.colorcet('C7'), colmaps.colorcet('C1')};
depthLimits = [-300 850];

features = {'direction', 'orientation'};
featureNames = {'dir', 'ori'};
[~, ~, mouseID] = unique({tuningData.animal});

%% Mean preference histograms across all datasets
% determine histograms for each dataset
dirHists = NaN(length(dirBins), length(tuningData));
oriHists = NaN(length(oriBins), length(tuningData));
for session = 1:length(tuningData)
    if sum(tuningData(session).dirTuned) >= minUnits
        pref = tuningData(session).dirPreferences(tuningData(session).dirTuned);
        dirHists(:,session) = histcounts(pref, dirEdges);
    end
    if sum(tuningData(session).oriTuned)
        pref = tuningData(session).oriPreferences(tuningData(session).oriTuned);
        oriHists(:,session) = histcounts(pref, oriEdges);
    end
end
indDir = ~any(isnan(dirHists), 1);
indOri = ~any(isnan(oriHists), 1);
% merge counts in first and last bars (around 0 and 360 or 180 deg)
dirHists(1,:) = dirHists(1,:) + dirHists(end,:);
dirHists(end,:) = dirHists(1,:);
oriHists(1,:) = oriHists(1,:) + oriHists(end,:);
oriHists(end,:) = oriHists(1,:);
% normalize histograms (to sum 1, without last bin)
dirHists = dirHists ./ sum(dirHists(1:end-1,:),1);
oriHists = oriHists ./ sum(oriHists(1:end-1,:),1);
% interpolate/smooth histograms
dirHistsSmooth = interp1([-30 dirBins 390], ...
    dirHists([end-1 1:end 2], indDir), dirBinsFine, 'pchip');
oriHistsSmooth = interp1([-15 oriBins 195], ...
    oriHists([end-1 1:end 2], indOri), oriBinsFine, 'pchip');

% plot direction
figure('Position', glob.figPositionDefault)
hold on
m = mean(dirHists, 2, "omitnan");
mSm = mean(dirHistsSmooth, 2, "omitnan");
sSm = std(dirHistsSmooth, 0, 2, "omitnan") ./ sqrt(sum(indDir));
fill(dirBinsFine([1:end end:-1:1]), [mSm-sSm; flip(mSm+sSm)], 'k', ...
    "FaceColor", 'k', "FaceAlpha", 0.5, "EdgeColor", "none")
plot(dirBinsFine, mSm, 'Color', 'k', "LineWidth", 1);
plot(dirBins, m, '.', 'Color', 'k', "MarkerSize", 30)
set(gca, "Box", "off", "XTick", 0:90:360)
xlim([-10 370])
ylim([0 0.2])
xlabel('Direction (deg)')
ylabel('Proportion of units')
title(sprintf('Direction tuning (%d sessions)', sum(indDir)))
io.saveFigure(gcf, fPlots, 'tuning_directionPrefHist');

% plot orientation
figure('Position', glob.figPositionDefault)
hold on
m = mean(oriHists, 2, "omitnan");
mSm = mean(oriHistsSmooth, 2, "omitnan");
sSm = std(oriHistsSmooth, 0, 2, "omitnan") ./ sqrt(sum(indOri));
fill(oriBinsFine([1:end end:-1:1]), [mSm-sSm; flip(mSm+sSm)], 'k', ...
    "FaceColor", 'k', "FaceAlpha", 0.5, "EdgeColor", "none")
plot(oriBinsFine, mSm, 'Color', 'k', "LineWidth", 1);
plot(oriBins, m, '.', 'Color', 'k', "MarkerSize", 30)
set(gca, "Box", "off", "XTick", 0:45:180)
xlim([-10 190])
ylim([0 0.2])
xlabel('Orientation (deg)')
ylabel('Proportion of units')
title(sprintf('Orientation tuning (%d sessions)', sum(indOri)))
io.saveFigure(gcf, fPlots, 'tuning_orientationPrefHist');

%% Preferences across depth (separately for each recording)
head = {'on', 'off'};
for feat = 1:2 % direction and orientation preferences
    figure('Position', [100 100 1090 650])
    hold on
    for session = 1:length(tuningData)
        units = find(tuningData(session).([featureNames{feat} 'Tuned']));
        X = ones(length(units),1) .* (session * 3 * scale);
        Y = tuningData(session).depth(units) - ...
            tuningData(session).SO_depth;
        angles = tuningData(session).([featureNames{feat} 'Preferences'])(units);
        if feat == 1 % direction
            [U, V] = pol2cart(deg2rad(angles), ...
                ones(length(units),1) .* scale);
        else % orientation
            % turn angles by 90 degrees to orientation of grating, rather
            % than direction of movement
            angles_plot = mod(angles + 90, 180);
            [U, V] = pol2cart(deg2rad(angles_plot), ...
                ones(length(units),1) .* scale);
        end
        X = X - 0.5 .* U;
        Y = Y - 0.5 .* V;

        c = colors{feat}(floor(angles ./ (360 / feat) .* 256) + 1, :);
        for unit = 1:length(units)
            quiver(X(unit), Y(unit), U(unit), V(unit), "off", ...
                'Color', c(unit,:), 'LineWidth', 2, ...
                'ShowArrowHead', head{feat}, 'MaxHeadSize', 8);
        end
    end
    plot([1.5 (length(tuningData) * 3 + 1.5)] .* scale, [0 0], 'k')
    axis image
    set(gca, "YDir", "reverse", "XTick", ...
        (1:length(tuningData)) .* (3 * scale), "XTickLabel", mouseID)
    ylim(depthLimits)
    xlabel('Animal ID')
    ylabel('Depth from SGS-SO border (in um)')
    colormap(colors{feat})
    if feat == 1 % direction
        colorbar("Ticks", (0:30:360) ./ 360, ...
            "TickLabels", {'F','','','D','','','B','','','U','','','F'})
    else % orientation
        colorbar("Ticks", (0:15:180) ./ 180, ...
            "TickLabels", {'V','','','','','','H','','','','','','V'})
    end
    title(sprintf('Preferred %ss', features{feat}))
    io.saveFigure(gcf, fPlots, ...
        sprintf('tuning_preference_%sAcrossSCDepth', features{feat}))
end

%% Selectivity across depth (pooled across recordings)
for feat = 1:2 % direction and orientation selectivities
    selectivities = [];
    depths_relative = [];
    for session = 1:length(tuningData)
        units = tuningData(session).([featureNames{feat} 'Tuned']);
        selectivities = [selectivities; ...
            tuningData(session).([featureNames{feat} 'Sel'])(units)];
        depths_relative = [depths_relative; ...
            tuningData(session).depth(units) - ...
            tuningData(session).SO_depth];
    end
    figure('Position', [100 100 300 650])
    hold on
    scatter(selectivities, depths_relative, 20, [1 1 1] .* 0.7, 'filled')
    [d, order] = sort(depths_relative);
    s = selectivities(order);
    f = fit(double(d), s, 'smoothingspline', 'SmoothingParam', 0.0001);
    xd = linspace(d(1), d(end), 100);
    plot(f(xd), xd, 'Color', [1 1 1] .* 0.3, 'LineWidth', 2)
    plot([0 1], [0 0], 'k')
    set(gca, "YDir", "reverse")
    xlabel(sprintf('%s selectivity', features{feat}))
    ylabel('Depth from SGS-SO border (in um)')
    title(sprintf('%s selectivity', features{feat}))
    xlim([0 1])
    ylim(depthLimits)
    io.saveFigure(gcf, fPlots, ...
        sprintf('tuning_selectivity_%sAcrossSCDepth', features{feat}))
end