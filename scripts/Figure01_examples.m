function Figure01_examples(folders, sets, ex)
% Example mean frames, calcium traces and tuning curves
buffer = 1; % in sec (before and after stim period)
for s = 1:2 % boutons and neurons
    str = sets{s};
    f = fullfile(folders.data, str, ex{s,1}, ex{s,2});
    % load data
    calc = io.getCalciumData(f);
    stim = io.getGratingInfo(f, 'gratingsDrifting');
    krnlFits = io.getStimResponseFits(f, 'gratingsDrifting');
    [dirTuning, oriTuning] = io.getTuningResults(f, 'gratingsDrifting');
    recInfo = io.getRecordingInfo(f);
    % stim. duration and non-gray stimuli
    stimDur = median(diff(stim.times,1,2));
    validStims = find(~isnan(stim.directions));

    units = ex{s,4};
    % align calcium trace to stimuli
    [alignedTrace, t] = traces.getAlignedTraces( ...
        calc.traces(:,units), calc.time, ...
        stim.times(:,1), [-buffer stimDur+buffer]);
    % align predicted trace to stimuli
    alignedPredicted = traces.getAlignedTraces( ...
        krnlFits.prediction(:,units), krnlFits.time_prediction, ...
        stim.times(:,1), [-buffer stimDur+buffer]);
    % loop over all examples
    for iUnit = 1:length(units)
        mini = min(alignedTrace(:,:,iUnit), [], "all");
        maxi = max(alignedTrace(:,:,iUnit), [], "all");
        f1 = figure('Position',[3 570 1915 420]);
        tiledlayout(1, length(validStims), "TileSpacing", "tight", ...
            "Padding", "tight")
        % loop over all stimuli
        for st = 1:length(validStims)
            indSt = stim.ids == validStims(st);
            nexttile
            hold on
            % mark time of stimulus
            fill([0 0 stimDur stimDur], [mini maxi maxi mini], 'k', ...
                'FaceColor', [1 1 1].*0.9, 'EdgeColor', 'none')
            % single trial traces
            plot(t, alignedTrace(:,indSt,iUnit), ...
                "Color", [1 1 1].*0.5);
            % mean predicted traces (from kernel fit)
            plot(t, mean(alignedPredicted(:,indSt,iUnit), 2, "omitnan"), ...
                'k', "LineWidth", 1)
            set(gca, "Box", "off")
            xlim(t([1 end]))
            ylim([mini maxi])
            if st == 1
                xlabel('Time (s)')
                ylabel('\DeltaF/F')
            end
            title(sprintf('%d deg', stim.directions(validStims(st))))
        end
        sgtitle(sprintf('ROI %03d: DS = %.2f (p=%.3f), OS = %.2f (p=%.3f)', ...
            units(iUnit), dirTuning.selectivity(units(iUnit)), ...
            dirTuning.pValue(units(iUnit)), oriTuning.selectivity(units(iUnit)), ...
            oriTuning.pValue(units(iUnit))))

        % plot fitted kernel
        mini = min(krnlFits.kernel(:,units(iUnit)));
        f2 = figure('Position', [150 150 475 325]);
        hold on
        fill([0 0 stimDur stimDur], [mini 1 1 mini], 'k', ...
                'FaceColor', [1 1 1].*0.9, 'EdgeColor', 'none')
        plot(krnlFits.time_kernel, krnlFits.kernel(:,units(iUnit)), ...
            'k', "LineWidth", 1)
        set(gca, "Box", "off")
        xlim(krnlFits.time_kernel([1 end]))
        ylim([mini 1])
        xlabel('Time (s)')
        title(sprintf('ROI %03d: kernel', units(iUnit)))

        % plot direction tuning curve
        amps = krnlFits.amplitudes(:, validStims([1:end 1]), units(iUnit));
        drct = [stim.directions(validStims); 360]';
        mini = min(amps, [], "all");
        maxi = max(amps, [], "all");
        range = maxi - mini;
        f3 = figure('Position', [700 150 475 325]);
        hold on
        plot(drct + randn(size(amps)).*3, ...
            amps, '.', 'Color', [1 1 1].*0.5, 'MarkerSize', 10);
        plot(drct, mean(amps,1,'omitnan'), ...
            'k.-', 'MarkerSize', 30, 'LineWidth', 1);
        plot(dirTuning.preference(units(iUnit)), maxi, 'v', 'MarkerSize', 8, ...
            'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'none');
        set(gca, 'box', 'off', 'XTick', drct(1:3:end))
        xlim([-10 370])
        ylim([min([0 mini-0.05*range]) maxi+0.05*range])
        xlabel('Direction (deg)')
        ylabel('\DeltaF/F (kernel amplitude)')
        title(sprintf('ROI %03d: Direction tuning', units(iUnit)))

        % plot orientation tuning curve
        amps = krnlFits.amplitudes(:, validStims, units(iUnit));
        amps = reshape(amps, size(amps,1), size(amps,2)/2, 2);
        amps = reshape(permute(amps, [1 3 2]), size(amps,1)*2, []);
        amps = amps(:,[1:end 1]);
        ortn = [stim.directions(validStims(1:length(validStims)/2)); 180]';
        mini = min(amps, [], "all");
        maxi = max(amps, [], "all");
        range = maxi - mini;
        f4 = figure('Position', [1250 150 475 325]);
        hold on
        plot(ortn + randn(size(amps)).*3, ...
            amps, '.', 'Color', [1 1 1].*0.5, 'MarkerSize', 10);
        plot(ortn, mean(amps,1,'omitnan'), ...
            'k.-', 'MarkerSize', 30, 'LineWidth', 1);
        plot(oriTuning.preference(units(iUnit)), maxi, 'v', 'MarkerSize', 8, ...
            'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'none');
        set(gca, 'box', 'off', 'XTick', ortn(1:3:end))
        xlim([-10 190])
        ylim([min([0 mini-0.05*range]) maxi+0.05*range])
        xlabel('Orientation (deg)')
        ylabel('\DeltaF/F (kernel amplitude)')
        title(sprintf('ROI %03d: Orientation tuning', units(iUnit)))

        io.saveFigure(f1, fPlot, sprintf('example_%s_stimTraces_%s_%s_%03d', ...
            str, ex{s,1}, ex{s,2}, units(iUnit)))
        io.saveFigure(f2, fPlot, sprintf('example_%s_kernel_%s_%s_%03d', ...
            str, ex{s,1}, ex{s,2}, units(iUnit)))
        io.saveFigure(f3, fPlot, sprintf('example_%s_dirTuningCurve_%s_%s_%03d', ...
            str, ex{s,1}, ex{s,2}, units(iUnit)))
        io.saveFigure(f4, fPlot, sprintf('example_%s_oriTuningCurve_%s_%s_%03d', ...
            str, ex{s,1}, ex{s,2}, units(iUnit)))
    end

    % plot ROI masks with ROI IDs
    masks = NaN(size(recInfo.roiMasks));
    masks(units,:) = recInfo.roiMasks(units,:);
    b = recInfo.fovBoundaries(ex{s,3},:);
    map = spatial.getROIMaskImage(masks, recInfo.fovPix(ex{s,3},:), b);
    colors = spatial.plotROIMaskImage(map, masks, true);
    io.saveFigure(gcf, fPlot, sprintf('example_%s_roiMasks_%s_%s_plane%02d', ...
        str, ex{s,1}, ex{s,2}, ex{s,3}))

    % plot ROI masks on mean image
    im = squeeze(recInfo.meanFrame(ex{s,3}, b(1):b(2), b(3):b(4)));
    imMasks = spatial.mergeImageWithMasks(im, map, colors);
    figure
    imshow(imMasks)
    set(gcf, 'Position', [680 50 1050 945])
    io.saveFigure(gcf, fPlot, sprintf('example_%s_roiMasksOnImage_%s_%s_plane%02d', ...
        str, ex{s,1}, ex{s,2}, ex{s,3}))

    % plot mean frame
    im = (im - min(im,[],"all"));
    im = im ./ max(im,[],"all");
    figure('Position', [680 50 1050 945])
    imagesc(imadjust(im))
    colormap(colmaps.getGCaMPMap)
    axis image off
    io.saveFigure(gcf, fPlot, sprintf('example_%s_meanFrame_%s_%s_plane%02d', ...
        str, ex{s,1}, ex{s,2}, ex{s,3}))
end