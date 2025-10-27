function Figure06_plotTuningData(tuningData, fPlots, glob)
scale = 30;
colors = {colmaps.colorcet('C7'), colmaps.colorcet('C1')};

features = {'direction', 'orientation'};
featureNames = {'dir', 'ori'};
[~, ~, mouseID] = unique({tuningData.animal});

% Preferences across depth (separately for each recording)
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

% Selectivity across depth (pooled across recordings)
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
    f = fit(d, s, 'smoothingspline', 'SmoothingParam', 0.0001);
    xd = linspace(d(1), d(end), 100);
    plot(f(xd), xd, 'Color', [1 1 1] .* 0.3, 'LineWidth', 2)
    plot([0 1], [0 0], 'k')
    set(gca, "YDir", "reverse")
    xlabel(sprintf('%s selectivity', features{feat}))
    ylabel('Depth from SGS-SO border (in um)')
    title(sprintf('%s selectivity', features{feat}))
    xlim([0 1])
    ylim(d([1 end]))
    io.saveFigure(gcf, fPlots, ...
        sprintf('tuning_selectivity_%sAcrossSCDepth', features{feat}))
end

%% Plot direction vs orientation preference and selectivity
tuned = [cat(1, tuningData.dirTuned), cat(1, tuningData.oriTuned)];
dirPreferences = cat(1, tuningData.dirPreferences);
oriPreferences = cat(1, tuningData.oriPreferences);
dirSel = cat(1, tuningData.dirSel);
oriSel = cat(1, tuningData.oriSel);

% Direction vs orientation preference scatterplot
figure('Position', glob.figPositionDefault)
hold on
plot([0 180], [0 180], 'Color', [1 1 1].*0.5)
plot([180 360], [0 180], 'Color', [1 1 1].*0.5)
scatter(dirPreferences(all(tuned,2)), oriPreferences(all(tuned,2)), 15, ...
    'k', 'filled')
axis equal
xlim([-10 370])
ylim([-10 190])
set(gca, "Box", "off", "XTick", 0:90:360, "YTick", 0:90:180)
xlabel('Direction (deg)')
ylabel('Orientation (deg)')
title(sprintf('n = %d', sum(all(tuned,2))))
io.saveFigure(gcf, fPlots, 'tuning_preference_dirVsOriScatter');

% DS vs OS scatterplot
figure('Position', glob.figPositionDefault)
h = gscatter(dirSel(any(tuned,2)), oriSel(any(tuned,2)), ...
    tuned(any(tuned,2),1) + 2*tuned(any(tuned,2),2), ...
    repmat([0.4 0.8 0]', 1, 3), [], 15);
% hold on
% scatter(dirSel(indExamples), oriSel(indExamples), 40, ...
%     lines(length(indExamples)), "filled")
l = legend(h, 'DS', 'OS', 'DS & OS', "Location", "bestoutside");
l.Box = "off";
axis padded equal
mini = -0.05;
maxi = 1.05;
axis([mini maxi mini maxi])
set(gca, "Box", "off")
xlabel('Direction selectivity')
ylabel('Orientation selectivity')
io.saveFigure(gcf, fPlots, 'tuning_selectivity_dirVsOriScatter');