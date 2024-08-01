%% Folders
getFolders;

%% Parameters
sets = {'boutons', 'neurons'};
maxP = 0.05; % p-value threshold for direction/orientation selectivity

%% Examples
ex = cell(2,4); % rows: (1) bouton, (2) neuron
% good boutons
% ex(1,:) = {'SS078', '2017-09-28', 1, [1 3 4 5 6 8 9 10 12 14 15 16 17 18 19 20 21 23 24 26 27 29 31 32 33 34 35 36 38 39 40 42 45 46 48 49 50 51 52 57 58 60 62 63 64 71 72 74 80 83 95 99 100 105 107 109 110 123 125 126 128 135 142 145 146 150 165 168 175 176 180 181 186 187 192 196 200 201 205 207 208 217 224 229 239 242 252 257 262 269 276]};
% final selection
% ex(1,:) = {'SS078', '2017-09-28', 1, [40 46 48  51 62 192  42  74 205]};
ex(1,:) = {'SS078', '2017-09-28', 1, [40 48  51 62 192  42]};
% best OS boutons
% ex(1,:) = {'SS078', '2017-09-28', 1, [18 45 150]};
% other good datasets
% ex(1,:) = {'SS078', '2017-10-05', 1, []};
% ex(1,:) = {'SS077', '2017-10-03', 1, []};
% good neurons
% ex(2,:) = {'SS044', '2015-04-28', 3, [227 231 236 240 245 252 256 258 268 281 285 293 313 321 328 369 378 379 380 383 385 393 408 425 426 430]};
% best neurons (mix of DS and OS)
% ex(2,:) = {'SS044', '2015-04-28', 3, [231 240 245 252 256 258 268 281 293 313 321 328 369 378 379 380 385 393]};
% best DS neurons
% ex(2,:) = {'SS044', '2015-04-28', 3, [227 236 240 285 383 408 425 426 430]};
% best OS neurons
% ex(2,:) = {'SS044', '2015-04-28', 3, [227 236 240 285 383 408 425 426 430]};
% final selection (2 triples of neighbours)
ex(2,:) = {'SS044', '2015-04-28', 3, [328 378 369 236 245 258]};

%% Add paths
addpath(genpath(fullfile(folders.tools, 'npy-matlab')))
addpath(fullfile(folders.repo))

%% For all plots
fPlot = fullfile(folders.plots, 'Figure01');
if ~isfolder(fPlot)
    mkdir(fPlot)
end

%% Example mean frames, calcium traces and tuning curves
buffer = 1; % in sec (before and after stim period)
for s = 1:2
    str = sets{s};
    f = fullfile(folders.data, str, ex{s,1}, ex{s,2});
    calc = io.getCalciumData(f);
    stim = io.getGratingInfo(f, 'gratingsDrifting');
    krnlFits = io.getStimResponseFits(f, 'gratingsDrifting');
    [dirTuning, oriTuning] = io.getTuningResults(f, 'gratingsDrifting');
    recInfo = io.getRecordingInfo(f);

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
        for st = 1:length(validStims)
            indSt = stim.ids == validStims(st);
            nexttile
            hold on
            fill([0 0 stimDur stimDur], [mini maxi maxi mini], 'k', ...
                'FaceColor', [1 1 1].*0.9, 'EdgeColor', 'none')
            plot(t, alignedTrace(:,indSt,iUnit), ...
                "Color", [1 1 1].*0.5);
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

        amps = krnlFits.amplitudes(:, validStims([1:end 1]), units(iUnit));
        drct = [stim.directions(validStims); 360]';
        mini = min(amps, [], "all");
        maxi = max(amps, [], "all");
        range = maxi - mini;
        f3 = figure('Position', [800 150 475 325]);
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

        io.saveFigure(f1, fPlot, sprintf('example_%s_stimTraces_%s_%s_%03d', ...
            str, ex{s,1}, ex{s,2}, units(iUnit)))
        io.saveFigure(f2, fPlot, sprintf('example_%s_kernel_%s_%s_%03d', ...
            str, ex{s,1}, ex{s,2}, units(iUnit)))
        io.saveFigure(f3, fPlot, sprintf('example_%s_tuningCurve_%s_%s_%03d', ...
            str, ex{s,1}, ex{s,2}, units(iUnit)))
    end

    % plot ROI masks with numbers
    masks = NaN(size(recInfo.roiMasks));
    masks(units,:) = recInfo.roiMasks(units,:);
    b = recInfo.fovBoundaries(ex{s,3},:);
    map = spatial.getROIMaskImage(masks, recInfo.fovPix(ex{s,3},:), b);
    colors = spatial.plotROIMaskImage(map, masks, true);
    io.saveFigure(gcf, fPlot, sprintf('example_%s_roiMasks_%s_%s_%03d', ...
        str, ex{s,1}, ex{s,2}, units(iUnit)))

    % plot ROI masks on mean image
    im = squeeze(recInfo.meanFrame(ex{s,3}, b(1):b(2), b(3):b(4)));
    imMasks = spatial.mergeImageWithMasks(im, map, colors);
    figure
    imshow(imMasks)
    set(gcf, 'Position', [680 50 1050 945])
    io.saveFigure(gcf, fPlot, sprintf('example_%s_roiMasksOnImage_%s_%s_%03d', ...
        str, ex{s,1}, ex{s,2}, units(iUnit)))
    % plot mean frame
    im = (im - min(im,[],"all"));
    im = im ./ max(im,[],"all");
    figure('Position', [680 50 1050 945])
    imagesc(imadjust(im))
    colormap(colmaps.getGCaMPMap)
    axis image off
    io.saveFigure(gcf, fPlot, sprintf('example_%s_meanFrame_%s_%s_%03d', ...
        str, ex{s,1}, ex{s,2}, units(iUnit)))
end

%% Population direction tuning curves
for s = 1:2 % neurons and boutons
    subjDirs = dir(fullfile(folders.data, sets{s}, 'SS*'));
    oriCurves = [];
    dirCurves = [];
    stimAll = 0:30:330;
    for subj = 1:length(subjDirs) % animals
        name = subjDirs(subj).name;
        dateDirs = dir(fullfile(folders.data, sets{s}, name, '2*'));
        for dt = 1:length(dateDirs) %dates
            date = dateDirs(dt).name;
            f = fullfile(folders.data, sets{s}, name, date);

            stim = io.getGratingInfo(f, 'gratingsDrifting');
            krnlFits = io.getStimResponseFits(f, 'gratingsDrifting');
            [dirTuning, oriTuning] = io.getTuningResults(f, 'gratingsDrifting');

        end
    end
end

%% Population orientation tuning curves

%% Total direction preferences histogram

%% Total orientation preferences histogram

%% Direction versus orientation selectivity indices