% Population tuning curves, preference histograms, DS vs OS
dirBinsCoarse = 0:30:360;
dirEdges = (0:30:390) - 15;
dirBinsFine = 0:5:360;
oriBinsCoarse = 0:30:180;
oriEdges = (0:30:210) - 15;
oriBinsFine = 0:5:180;

for s = 1:2 % boutons and neurons
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
            % load data
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

    % normalize tuning curves to range [0 1]
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
    dirCurves = interp1(dirBinsCoarse, dirCurves', dirBinsFine, "pchip")';
    oriCurves = interp1(oriBinsCoarse, oriCurves', oriBinsFine, "pchip")';

    % plot direction tuning curves
    figure('Position', [30 50 485 945])
    imagesc(dirCurves, [0 1])
    colormap(flip(gray, 1))
    set(gca, "Box", "off", "XTick", 1:6:length(dirBinsFine), "XTickLabel", dirBinsCoarse)
    xlabel('Direction (deg)')
    ylabel(sets{s})
    io.saveFigure(gcf, fPlot, sprintf('tuning_%s_directionHeatmap', sets{s}));
    % plot orientation tuning curves
    figure('Position', [600 50 485 945])
    imagesc(oriCurves, [0 1])
    colormap(flip(gray, 1))
    set(gca, "Box", "off", "XTick", 1:6:length(oriBinsFine), "XTickLabel", oriBinsCoarse)
    xlabel('Orientation (deg)')
    ylabel(sets{s})
    io.saveFigure(gcf, fPlot, sprintf('tuning_%s_orientationHeatmap', sets{s}));

    % plot direction preference histogram
    figure
    n1 = histcounts(dirPreferences(R2>minR2 & dirTuned & oriTuned), dirEdges);
    n2 = histcounts(dirPreferences(R2>minR2 & dirTuned & ~oriTuned), dirEdges);
    b = bar(dirBinsCoarse, [n1' n2'], 'stacked');
    b(1).FaceColor = 'k';
    b(2).FaceColor = [.5 .5 .5];
    l = legend(b, sprintf('DS & OS (%d)', sum(R2>minR2 & dirTuned & oriTuned)), ...
        sprintf('DS (%d)', sum(R2>minR2 & dirTuned & ~oriTuned)));
    l.Box = "off";
    xlim([-20 380])
    set(gca, "Box", "off", "XTick", 0:90:360)
    xlabel('Direction (deg)')
    ylabel(sets{s})
    io.saveFigure(gcf, fPlot, sprintf('tuning_%s_directionPrefHist', sets{s}));
    % plot orientation preference histogram
    figure
    n1 = histcounts(oriPreferences(R2>minR2 & oriTuned & dirTuned), oriEdges);
    n2 = histcounts(oriPreferences(R2>minR2 & oriTuned & ~dirTuned), oriEdges);
    b = bar(oriBinsCoarse, [n1' n2'], 'stacked');
    b(1).FaceColor = 'k';
    b(2).FaceColor = [.5 .5 .5];
    l = legend(b, sprintf('OS & DS (%d)', sum(R2>minR2 & dirTuned & oriTuned)), ...
        sprintf('OS (%d)', sum(R2>minR2 & ~dirTuned & oriTuned)));
    l.Box = "off";
    xlim([-20 200])
    set(gca, "Box", "off", "XTick", 0:90:180)
    xlabel('Orientation (deg)')
    ylabel(sets{s})
    io.saveFigure(gcf, fPlot, sprintf('tuning_%s_orientationPrefHist', sets{s}));

    % plot direction preference histogram per dataset
    dirHists = NaN(length(dirBinsCoarse), max(dataset));
    oriHists = NaN(length(oriBinsCoarse), max(dataset));
    nDir = histcounts(dataset(R2>minR2 & dirTuned), 0.5:max(dataset)+1);
    nOri = histcounts(dataset(R2>minR2 & oriTuned), 0.5:max(dataset)+1);
    for c = 1:max(dataset)
        pref = dirPreferences(R2>minR2 & dirTuned & c==dataset);
        dirHists(:,c) = histcounts(pref, dirEdges);
        pref = oriPreferences(R2>minR2 & oriTuned & c==dataset);
        oriHists(:,c) = histcounts(pref, oriEdges);
    end
    % interpolate/smooth and normalize histograms (to sum 1)
    dirHists = interp1(dirBinsCoarse, dirHists, ...
        dirBinsFine, 'pchip');
    dirHists = dirHists ./ sum(dirHists,1);
    oriHists = interp1(oriBinsCoarse, oriHists, ...
        oriBinsFine, 'pchip');
    oriHists = oriHists ./ sum(oriHists,1);
    figure
    plot(dirBinsFine, dirHists, "LineWidth", 2)
    l = legend(num2str(nDir'), "Location", "bestoutside");
    l.Title.String = sprintf('#%s', sets{s});
    l.Box = "off";
    xlim([-20 380])
    set(gca, "Box", "off", "ColorOrder", turbo(max(dataset)), ...
        "XTick", 0:90:360)
    xlabel('Direction (deg)')
    ylabel(sprintf('Proportion %s', sets{s}))
    io.saveFigure(gcf, fPlot, ...
        sprintf('tuning_%s_directionPrefHistPerDataset', sets{s}));
    % plot orientation preference histogram per dataset
    figure
    plot(oriBinsFine, oriHists, "LineWidth", 2)
    l = legend(num2str(nOri'), "Location", "bestoutside");
    l.Title.String = sprintf('#%s', sets{s});
    l.Box = "off";
    xlim([-20 200])
    set(gca, "Box", "off", "ColorOrder", turbo(max(dataset)), ...
        "XTick", 0:90:180)
    xlabel('Orientation (deg)')
    ylabel(sprintf('Proportion %s', sets{s}))
    io.saveFigure(gcf, fPlot, ...
        sprintf('tuning_%s_orientationPrefHistPerDataset', sets{s}));

    % plot DS vs OS scatterplot
    ind = R2>minR2 & (dirTuned | oriTuned);
    figure
    h = gscatter(dirSel(ind), oriSel(ind), dirTuned(ind) + ...
        2*oriTuned(ind), [], [], 15);
    l = legend(h, 'DS', 'OS', 'DS & OS', "Location", "bestoutside");
    l.Box = "off";
    axis padded equal
    mini = min(axis);
    maxi = max(axis);
    axis([mini maxi mini maxi])
    set(gca, "Box", "off")
    xlabel('Direction selectivity')
    ylabel('Orientation selectivity')
    title(sprintf('%s', sets{s}))
    io.saveFigure(gcf, fPlot, ...
        sprintf('tuning_%s_selectivityDirVsOri', sets{s}));
end