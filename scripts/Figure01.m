%% Folders
getFolders;

%% Parameters
sets = {'boutons', 'neurons'};
maxP = 0.05; % p-value threshold for response kernel and 
             % direction/orientation selectivity
minR2 = 0.02; % threshold for explained variance for response kernel

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
    curves = [];
    R2 = [];
    dirPreferences = [];
    oriPreferences = [];
    dirSel = [];
    oriSel = [];
    dirTuned = [];
    oriTuned = [];
    dataset = [];
    stimAll = 0:30:330;
    count = 1;
    for subj = 1:length(subjDirs) % animals
        name = subjDirs(subj).name;
        dateDirs = dir(fullfile(folders.data, sets{s}, name, '2*'));
        for d = 1:length(dateDirs) %dates
            date = dateDirs(d).name;
            f = fullfile(folders.data, sets{s}, name, date);
            % ignore session if stimulus was not presented
            if ~isfile(fullfile(f, '_ss_gratingsDrifting.intervals.npy'))
                continue
            end

            stim = io.getGratingInfo(f, 'gratingsDrifting');
            krnlFits = io.getStimResponseFits(f, 'gratingsDrifting');
            [dirTuning, oriTuning] = io.getTuningResults(f, 'gratingsDrifting');

            % check whether the typical stimulus directions were used in
            % experiment
            if ~all(ismember(stimAll, stim.directions))
                fprintf('WARNING: Grating directions are different from usual in %s %s!\n', ...
                    name, date)
            end
            % find indices of relevant stimuli
            indStims = find(ismember(stim.directions, stimAll));
            % find units that are responsive
            unitsResponsive = krnlFits.pValue < maxP;
            % determine tuning curves
            c = NaN(sum(unitsResponsive), length(stimAll));
            c(:,indStims) = squeeze(mean( ...
                krnlFits.amplitudes(:,indStims,unitsResponsive), 1, "omitnan"))';
            % append tuning curves of current dataset
            curves = [curves; c];
            R2 = [R2; krnlFits.R2(unitsResponsive)];
            dirPreferences = [dirPreferences; dirTuning.preference(unitsResponsive)];
            oriPreferences = [oriPreferences; oriTuning.preference(unitsResponsive)];
            dirSel = [dirSel; dirTuning.selectivity(unitsResponsive)];
            oriSel = [oriSel; oriTuning.selectivity(unitsResponsive)];
            dirTuned = [dirTuned; dirTuning.pValue(unitsResponsive) < maxP];
            oriTuned = [oriTuned; oriTuning.pValue(unitsResponsive) < maxP];
            dataset = [dataset; ones(sum(unitsResponsive),1) .* count];

            count = count + 1;
        end
    end
    % direction tuning curves: direction selective units (sorted by pref.
    % dir.) + gap + non-selective units
    dc = curves(dirTuned & R2>minR2,:);
    [~, order] = sort(dirPreferences(dirTuned & R2>minR2));
    dirCurves = [dc(order,:); zeros(100,length(stimAll))];
    dc = curves(~dirTuned & R2>minR2,:);
    [~, order] = sort(dirPreferences(~dirTuned & R2>minR2));
    dirCurves = [dirCurves; dc(order,:)];
    % determine orientation tuning curves
    oCurves = permute(curves, [2 3 1]);
    oCurves = reshape(oCurves, length(stimAll)/2, 2, []);
    oCurves = mean(oCurves,2);
    oCurves = permute(oCurves, [3 1 2]);
    % orientation tuning curves: orientation selective units (sorted by pref.
    % dir.) + gap + non-selective units
    oc = oCurves(oriTuned & R2>minR2,:);
    [~, order] = sort(oriPreferences(oriTuned & R2>minR2));
    oriCurves = [oc(order,:); NaN(100, size(oc,2))];
    oc = oCurves(~oriTuned & R2>minR2,:);
    [~, order] = sort(oriPreferences(~oriTuned & R2>minR2));
    oriCurves = [oriCurves; oc(order,:)];
    % normalize tuning curves
    dirCurves = dirCurves - min(dirCurves, [], 2);
    dirCurves = dirCurves ./ max(dirCurves, [], 2);
    dirCurves(isnan(dirCurves)) = 0;
    oriCurves = oriCurves - min(oriCurves, [], 2);
    oriCurves = oriCurves ./ max(oriCurves, [], 2);
    oriCurves(isnan(oriCurves)) = 0;
    % duplicate responses to 0 deg (for 360 deg and 180 deg)
    dirCurves = dirCurves(:,[1:end 1]);
    oriCurves = oriCurves(:,[1:end 1]);
    % upsample and smooth curves
    dirs1 = 0:30:360;
    dirs2 = 0:5:360;
    oris1 = 0:30:180;
    oris2 = 0:5:180;
    dirCurves = interp1(dirs1, dirCurves', dirs2, "pchip")';
    oriCurves = interp1(oris1, oriCurves', oris2, "pchip")';

    % plot direction tuning curves
    figure('Position', [150 50 560 785])
    imagesc(dirCurves, [0 1])
    colormap(flip(gray, 1))
    set(gca, "Box", "off", "XTick", 1:6:length(dirs2), "XTickLabel", dirs1)
    xlabel('Direction (deg)')
    ylabel(sets{s})
    io.saveFigure(gcf, fPlot, sprintf('tuning_%s_directionHeatmap', sets{s}));
    % plot orientation tuning curves
    figure('Position', [150 50 560 785])
    imagesc(oriCurves, [0 1])
    colormap(flip(gray, 1))
    set(gca, "Box", "off", "XTick", 1:6:length(oris2), "XTickLabel", oris1)
    xlabel('Direction (deg)')
    ylabel(sets{s})
    io.saveFigure(gcf, fPlot, sprintf('tuning_%s_orientationHeatmap', sets{s}));

    % plot direction preference histogram
    figure
    dirBins = 0:30:360;
    dirEdges = (0:30:390) - 15;
    n1 = histcounts(dirPreferences(R2>minR2 & dirTuned & oriTuned), dirEdges);
    n2 = histcounts(dirPreferences(R2>minR2 & dirTuned & ~oriTuned), dirEdges);
    b = bar(dirBins, [n1' n2'], 'stacked');
    b(1).FaceColor = 'k';
    b(2).FaceColor = [.5 .5 .5];
    xlim([-20 380])
    set(gca, "Box", "off", "XTick", 0:90:360)
    xlabel('Direction (deg)')
    ylabel(sets{s})
    legend(b, sprintf('DS & OS (%d)', sum(R2>minR2 & dirTuned & oriTuned)), ...
        sprintf('DS (%d)', sum(R2>minR2 & dirTuned & ~oriTuned)))
    io.saveFigure(gcf, fPlot, sprintf('tuning_%s_directionPrefHist', sets{s}));
    % plot orientation preference histogram
    figure
    oriBins = 0:30:180;
    oriEdges = (0:30:210) - 15;
    n1 = histcounts(oriPreferences(R2>minR2 & oriTuned & dirTuned), oriEdges);
    n2 = histcounts(oriPreferences(R2>minR2 & oriTuned & ~dirTuned), oriEdges);
    b = bar(oriBins, [n1' n2'], 'stacked');
    b(1).FaceColor = 'k';
    b(2).FaceColor = [.5 .5 .5];
    xlim([-20 200])
    set(gca, "Box", "off", "XTick", 0:90:180)
    xlabel('Orientation (deg)')
    ylabel(sets{s})
    legend(b, sprintf('OS & DS (%d)', sum(R2>minR2 & dirTuned & oriTuned)), ...
        sprintf('OS (%d)', sum(R2>minR2 & ~dirTuned & oriTuned)))
    io.saveFigure(gcf, fPlot, sprintf('tuning_%s_orientationPrefHist', sets{s}));

    % % plot direction preference histogram per dataset
    % dirBins = 0:30:360;
    % dirEdges = (0:30:390) - 15;
    % oriBins = 0:30:180;
    % oriEdges = (0:30:210) - 15;
    % dirHists = NaN(length(dirBins), max(dataset));
    % oriHists = NaN(length(oriBins), max(dataset));
    % for c = 1:max(dataset)
    %     pref = dirPreferences(R2>minR2 & dirTuned & c==dataset);
    %     dirHists(:,c) = histcounts(pref, dirEdges);
    %     pref = oriPreferences(R2>minR2 & oriTuned & c==dataset);
    %     oriHists(:,c) = histcounts(pref, oriEdges);
    % end
    % dirHists = dirHists ./ sum(dirHists,2);
    % oriHists = oriHists ./ sum(oriHists,2);
    % figure
    % plot(dirBins, dirHists', "LineWidth", 2)
    % set(gca, "ColorOrder", hsv(max(dataset)))

    % plot DS vs OS scatterplot
    figure
    hold on
    plot()
end

%% Total direction preferences histogram

%% Total orientation preferences histogram

%% Direction versus orientation selectivity indices