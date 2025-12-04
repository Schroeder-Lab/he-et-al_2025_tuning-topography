function Figure03_examples(folders, sets, ex, fPlots)

%% Parameters
% plotting RFs
edges_rf = [-135 -90 42 -20];
RFtypes = {'ON', 'OFF', 'ON+OFF'};
colEx = lines(4);

% for evaluation of receptive fields (significance/goodness)
minEV = 0.01; % minimum explained variance to plot RF
minPeak = 5; % minimum peak of RF (compared to noise) to plot RF

%% Plots of example data sets and units
for s = 1:2 % boutons and neurons
    str = sets{s};
    f = fullfile(folders.data, str, ex{s,1}, ex{s,2});

    % load data
    recInfo = io.getRecordingInfo(f);
    rfData = io.getRFFits(f);
    edges = rfData.edges;
    gridW = diff(edges(1:2)) / size(rfData.maps,3);
    gridH = -diff(edges(3:4)) / size(rfData.maps,2);
    rfGaussPars = rfData.fitParameters;
    EVs = rfData.EV;
    peakNoiseRatio = rfData.peaks;

    units = ex{s,4};

    % plot ROI masks with ROI IDs
    masks = NaN(size(recInfo.roiMasks));
    masks(units,:) = recInfo.roiMasks(units,:);
    b = recInfo.fovBoundaries(ex{s,3},:);
    map = spatial.getROIMaskImage(masks, recInfo.fovPix(ex{s,3},:), b);
    colors = spatial.plotROIMaskImage(map, masks, true);
    io.saveFigure(gcf, fPlots, sprintf('example_%s_roiMasks_%s_%s_plane%02d', ...
        str, ex{s,1}, ex{s,2}, ex{s,3}))

    % plot ROI masks on mean image
    im = squeeze(recInfo.meanFrame(ex{s,3}, b(1):b(2), b(3):b(4)));
    imMasks = spatial.mergeImageWithMasks(im, map, colors);
    figure
    imshow(imMasks)
    set(gcf, 'Position', [680 50 1050 945])
    io.saveFigure(gcf, fPlots, sprintf('example_%s_roiMasksOnImage_%s_%s_plane%02d', ...
        str, ex{s,1}, ex{s,2}, ex{s,3}))

    % plot mean frame
    im = (im - min(im,[],"all"));
    im = im ./ max(im,[],"all");
    figure('Position', [680 50 1050 945])
    imagesc(imadjust(im))
    colormap(colmaps.getGCaMPMap)
    axis image off
    io.saveFigure(gcf, fPlots, sprintf('example_%s_meanFrame_%s_%s_plane%02d', ...
        str, ex{s,1}, ex{s,2}, ex{s,3}))

    % example RFs
    for iUnit = units
        % rf_tmp: [rows x columns x subfield]
        rf_tmp = squeeze(mean(rfData.maps(iUnit,:,:,:,:),4));
        rf.plotRF(rf_tmp, rfGaussPars(iUnit,:), ...
            rfData.bestSubFields(iUnit), edges, edges_rf, gridW, gridH)
        sgtitle(sprintf('ROI %d (EV: %.3f, peak/noise: %.1f, %s)', ...
            iUnit, EVs(iUnit), peakNoiseRatio(iUnit), ...
            RFtypes{rfData.bestSubFields(iUnit)}))

        io.saveFigure(gcf, fPlots, sprintf('example_%s_RF_%s_%s_%03d', ...
            str, ex{s,1}, ex{s,2}, iUnit))
    end

    % RF outlines of all units in example dataset
    rf.plotRFOutlines(rfGaussPars, EVs, peakNoiseRatio, minEV, minPeak, ...
        units, edges_rf)
    n = sum(EVs >= minEV & peakNoiseRatio >= minPeak);
    title(sprintf('%s %s (n = %d)', ex{s,1}, ex{s,2}, n))
    io.saveFigure(gcf, fPlots, sprintf('example_%s_RFoutlines_%s_%s', ...
        str, ex{s,1}, ex{s,2}))

    % map RF to brain position for all units in example dataset
    % load data
    brainPos = recInfo.roiPositions; % (x,y,z) in microns
    brainPos(:,3) = [];
    rfPos = rfGaussPars(:, [2 4]); % (azimuth, elevation) in visual degrees
    outliers = readNPY(fullfile(f, '_ss_rf.outliers.npy'));
    model = readNPY(fullfile(f, "_ss_rfRetinotopy.model.npy"));
    fit_rfPos = @(x) rf.brainToRFPos(x, model(1), model(2), ...
        model(3), model(4), model(5));

    % only consider significant RFs
    valid = EVs >= minEV & peakNoiseRatio >= minPeak;
    valid = valid & ~outliers;

    % for all units with a significant RF, plot RF position and
    % colour code (1) horizontal and (2) vertical brain position
    brainLimits = reshape([floor(min(brainPos(valid,:),[],1)); ...
        ceil(max(brainPos(valid,:),[],1))], [], 1)';
    brainRange = [diff(brainLimits(1:2)) diff(brainLimits(3:4))];
    visualLimits = reshape([floor(min(rfPos(valid,:),[],1)); ...
        ceil(max(rfPos(valid,:),[],1))], [], 1)';
    [x, y] = meshgrid(linspace(brainLimits(1)-brainRange(1), ...
        brainLimits(2)+brainRange(1), 60), ...
        linspace(brainLimits(3)-brainRange(2), ...
        brainLimits(4)+brainRange(2), 60));
    rfXY = fit_rfPos([x(:); y(:)]);
    rfX = reshape(rfXY(1:3600),60,60);
    rfY = reshape(rfXY(3601:end),60,60);
    plotLimits = brainLimits + 0.1.*[-brainRange(1), brainRange(1), ...
        -brainRange(2), brainRange(2)];
    contourLimits = fit_rfPos([plotLimits([1 1 2 2]), ...
        plotLimits([3 4 3 4])]);

    figure('Position', [80 255 1105 420])
    tiledlayout(1,2)

    nexttile
    hold on
    cLimits = [min([contourLimits(1:4) visualLimits(1)]), ...
        max([contourLimits(1:4) visualLimits(2)])];
    [~,c] = contourf(x, y, rfX, floor(cLimits(1)):floor(cLimits(2)));
    clim(cLimits)
    c.LineStyle = "none";
    scatter(brainPos(valid,1), brainPos(valid,2), [], ...
        rfPos(valid,1), "filled", 'MarkerEdgeColor', 'k');
    for k = 1:length(units)
        iUnit = units(k);
        plot(brainPos(iUnit,1), brainPos(iUnit,2), 'o', ...
            'MarkerEdgeColor', colEx(k,:), 'LineWidth', 2)
    end
    colormap turbo
    c = colorbar;
    c.Label.String = 'RF azimuth';
    axis image
    axis(plotLimits)
    set(gca, "Box", "off", "YDir", "reverse")
    xlabel('Brain ML (\mum)')
    ylabel('Brain AP (\mum)')
    title(sprintf('Brain position vs RF azimuth (n = %d)', sum(valid)))

    nexttile
    hold on
    cLimits = [min([contourLimits(5:end) visualLimits(3)]), ...
        max([contourLimits(5:end) visualLimits(4)])];
    [~,c] = contourf(x, y, rfY, floor(cLimits(1)):floor(cLimits(2)));
    clim(cLimits)
    c.LineStyle = "none";
    scatter(brainPos(valid,1), brainPos(valid,2), [], ...
        rfPos(valid,2), "filled", 'MarkerEdgeColor', 'k')
    for k = 1:length(units)
        iUnit = units(k);
        plot(brainPos(iUnit,1), brainPos(iUnit,2), 'o', ...
            'MarkerEdgeColor', colEx(k,:), 'LineWidth', 2)
    end
    colormap turbo
    c = colorbar;
    c.Label.String = 'RF elevation';
    axis image
    axis(plotLimits)
    set(gca, "Box", "off", "YDir", "reverse")
    xlabel('Brain ML (\mum)')
    ylabel('Brain AP (\mum)')
    title('Brain position vs RF elevation')
    sgtitle(sprintf('%s %s', ex{s,1}, ex{s,2}), 'FontWeight', 'bold')
    io.saveFigure(gcf, fPlots, sprintf('example_%s_RF-brain-position_%s_%s', ...
        str, ex{s,1}, ex{s,2}))
end