function Figure03_examples(folders, sets, ex, fPlots)

%% Parameters
% plotting RFs
[cm_ON, cm_OFF] = colmaps.getRFMaps;
cms = cat(3, cm_ON, cm_OFF);
ellipse_x = linspace(-pi, pi, 100);
edges_rf = [-135 -90 42 -20];
titles = {'ON field','OFF field'};
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
        rf_tmp(:,:,2) = -rf_tmp(:,:,2);
        mx = max(abs(rf_tmp),[],"all");
        pars = rfGaussPars(iUnit,:);
        figure('Position', [75 195 1470 475])
        for sf = 1:2
            subplot(1,2,sf)
            % STA
            imagesc([edges(1)+gridW/2 edges(2)-gridW/2], ...
                [edges(3)-gridH/2 edges(4)+gridH/2], ...
                rf_tmp(:,:,sf),[-mx mx])
            hold on
            if ismember(rfData.bestSubFields(iUnit), [sf 3])
                line = '-';
            else
                line = ':';
            end
            % ellipse at 1 STD (x and y), not rotated, not shifted
            x = pars(3) * cos(ellipse_x);
            y = pars(5) * sin(ellipse_x);
            % rotate and shift ellipse
            x_rot = pars(2) + ...
                x .* cos(pars(6)) - ...
                y .* sin(pars(6));
            y_rot = pars(4) + ...
                x .* sin(pars(6)) + ...
                y .* cos(pars(6));
            plot(x_rot, y_rot, ['k' line], 'LineWidth', 2)
            axis image
            set(gca, 'box', 'off', 'YDir', 'normal')
            xlim(edges_rf([1 2]))
            ylim(edges_rf([4 3]))
            colormap(gca, cms(:,:,sf))
            title(titles{sf})
            colorbar
        end
        sgtitle(sprintf('ROI %d (EV: %.3f, peak/noise: %.1f, %s)', ...
            iUnit, EVs(iUnit), peakNoiseRatio(iUnit), ...
            RFtypes{rfData.bestSubFields(iUnit)}))

        io.saveFigure(gcf, fPlots, sprintf('example_%s_RF_%s_%s_%03d', ...
            str, ex{s,1}, ex{s,2}, iUnit))
    end

    % RF outlines of all units in example dataset
    figure
    hold on
    h = zeros(size(rfGaussPars,1),1);
    for iUnit = 1:size(rfGaussPars,1)
        if EVs(iUnit) < minEV || peakNoiseRatio(iUnit) < minPeak
            continue
        end
        % ellipse at 2 STD (x and y), not rotated, not shifted
        x = rfGaussPars(iUnit,3) * cos(ellipse_x);
        y = rfGaussPars(iUnit,5) * sin(ellipse_x);
        % rotate and shift ellipse
        x_rot = rfGaussPars(iUnit,2) + ...
            x .* cos(rfGaussPars(iUnit,6)) - ...
            y .* sin(rfGaussPars(iUnit,6));
        y_rot = rfGaussPars(iUnit,4) + ...
            x .* sin(rfGaussPars(iUnit,6)) + ...
            y .* cos(rfGaussPars(iUnit,6));
        h(iUnit) = plot(x_rot, y_rot, 'k');
    end
    for k = 1:length(units)
        iUnit = units(k);
        set(h(iUnit),'Color',colEx(k,:))
        uistack(h(iUnit),'top')
        set(h(iUnit),'LineWidth',2)
    end
    legend(h(units))
    axis image
    axis(edges_rf([1 2 4 3]))
    xlabel('Azimuth (visual degrees)')
    ylabel('Elevation (visual degrees)')
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