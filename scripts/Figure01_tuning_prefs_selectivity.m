function Figure01_tuning_prefs_selectivity(folders, sets, maxP, minR2, ex, fPlots)
% Population tuning curves, preference histograms, DS vs OS
minUnits = 20;

dirBinsCoarse = 0:30:360;
dirEdges = (0:30:390) - 15;
dirBinsFine = 0:1:360;
oriBinsCoarse = 0:30:180;
oriBinsFine = 0:1:180;
oriEdges = (0:15:195) - 7.5;
oriBins = 0:15:180;

figDir = figure;
hDir = [0 0];
hold on
figOri = figure;
hOri = [0 0];
hold on
figDS = figure;
hold on
figOS = figure;
hold on
numDir = [0 0];
numOri = [0 0];
colSets = lines(2);

medianDS = cell(1,2);
medianOS = cell(1,2);
fprintf('In response to drifting gratings:\n')
for s = 1:2 % boutons and neurons
    fprintf('  %s:\n', sets{s})
    subjDirs = dir(fullfile(folders.data, sets{s}, 'SS*'));
    curves = [];
    R2 = [];
    dirPreferences = [];
    oriPreferences = [];
    dirSel = [];
    oriSel = [];
    dirTuned = false(0,0);
    oriTuned = false(0,0);
    depth = [];
    isGad = NaN(0,0);
    dataset = [];
    totalN = [];
    stimAll = 0:30:330;
    animals = 1;
    sessions = 1;
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
            recData = io.getRecordingInfo(f);
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
            unitsResponsive = krnlFits.pValue < maxP & krnlFits.R2 > minR2;
            % determine tuning curves
            c = NaN(sum(unitsResponsive), length(stimAll));
            c(:,indStims) = squeeze(mean( ...
                krnlFits.amplitudes(:,indStims,unitsResponsive), 1, "omitnan"))';
            % append tuning curves of current dataset
            curves = [curves; c .* dirTuning.responseSign(unitsResponsive)];
            R2 = [R2; krnlFits.R2(unitsResponsive)];
            dirPreferences = [dirPreferences; dirTuning.preference(unitsResponsive)];
            oriPreferences = [oriPreferences; oriTuning.preference(unitsResponsive)];
            dirSel = [dirSel; dirTuning.selectivity(unitsResponsive)];
            oriSel = [oriSel; oriTuning.selectivity(unitsResponsive)];
            dirTuned = [dirTuned; dirTuning.pValue(unitsResponsive) < maxP];
            oriTuned = [oriTuned; oriTuning.pValue(unitsResponsive) < maxP];
            depth = [depth; recData.surface - ...
                recData.roiPositions(unitsResponsive, 3)];
            isGad = [isGad; recData.isInhibitory(unitsResponsive)];
            dataset = [dataset; ones(sum(unitsResponsive),1) .* sessions];
            totalN = [totalN; length(krnlFits.pValue)];

            sessions = sessions + 1;

            if strcmp(name, ex{s,1}) && strcmp(date, ex{s,2})
                indExamples = NaN(length(ex{s,4}),1);
                unitsResponsive = find(unitsResponsive);
                n = length(dataset) - length(unitsResponsive);
                for k = 1:length(indExamples)
                    indExamples(k) = n + find(ex{s,4}(k) == unitsResponsive);
                end
            end
        end
        animals = animals + 1;
    end
    dirPreferences = mod(dirPreferences, 360);
    oriPreferences = mod(oriPreferences, 180);
    
    fprintf('    %d were recorded across %d sessions from %d mice\n', ...
        sum(totalN), sessions-1, animals-1)
    fprintf('    Explained variance of kernel-based predictions (mean +- STD): %.4f +- %.4f\n', ...
        mean(R2(R2 > minR2)), std(R2(R2 > minR2)))
    fprintf('    Of %d responsives, %d (%.1f%%) were tuned to direction, %d (%.1f%%) were tuned to orientation\n', ...
        sum(R2 > minR2), sum(dirTuned), sum(dirTuned)/length(dataset)*100, ...
        sum(oriTuned), sum(oriTuned)/length(dataset)*100)

    indEx = zeros(size(dataset));
    for k = 1:length(indExamples)
        indEx(indExamples(k)) = k;
    end

    %% Heatmaps of tuning curves
    % direction tuning curves: direction selective units (sorted by pref.
    % dir.) + gap + non-selective units
    dc = curves(dirTuned,:);
    [~, order] = sort(dirPreferences(dirTuned));
    dirCurves = [dc(order,:); zeros(100,length(stimAll))];
    indExDir = indEx(dirTuned);
    indExOrderedDir = NaN(size(indExamples));
    for k = 1:length(indExamples)
        indExOrderedDir(k) = find(order == find(indExDir == k));
    end
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
    oc = oCurves(oriTuned,:);
    [~, order] = sort(oriPreferences(oriTuned));
    oriCurves = [oc(order,:); NaN(100, size(oc,2))];
    indExOri = indEx(oriTuned);
    indExOrderedOri = NaN(size(indExamples));
    for k = 1:length(indExamples)
        indExOrderedOri(k) = find(order == find(indExOri == k));
    end
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

    % plot direction tuning curves (heatmap)
    figure('Position', [30 50 485 945])
    imagesc(dirCurves, [0 1])
    colormap(flip(gray, 1))
    hold on
    h = plot(3, indExOrderedDir, '>');
    legend(h, num2str(ex{s,4}'))
    set(gca, "Box", "off", "XTick", 1:6:length(dirBinsFine), "XTickLabel", dirBinsCoarse)
    xlabel('Direction (deg)')
    ylabel(sets{s})
    io.saveFigure(gcf, fPlots, sprintf('tuning_%s_directionHeatmap', sets{s}));

    % plot orientation tuning curves (heatmap)
    figure('Position', [600 50 485 945])
    imagesc(oriCurves, [0 1])
    colormap(flip(gray, 1))
    hold on
    h = plot(3, indExOrderedOri, '>');
    legend(h, num2str(ex{s,4}'))
    set(gca, "Box", "off", "XTick", 1:6:length(oriBinsFine), "XTickLabel", oriBinsCoarse)
    xlabel('Orientation (deg)')
    ylabel(sets{s})
    io.saveFigure(gcf, fPlots, sprintf('tuning_%s_orientationHeatmap', sets{s}));

    %% Preference histograms for all datasets
    % direction
    figure
    n1 = histcounts(dirPreferences(dirTuned & oriTuned), dirEdges);
    n2 = histcounts(dirPreferences(dirTuned & ~oriTuned), dirEdges);
    % merge counts in first and last bars (around 0 and 360 deg)
    n1(1) = n1(1) + n1(end);
    n1(end) = n1(1);
    n2(1) = n2(1) + n2(end);
    n2(end) = n2(1);
    b = bar(dirBinsCoarse, [n1' n2'], 'stacked');
    b(1).FaceColor = 'k';
    b(2).FaceColor = [.5 .5 .5];
    ylim([0 400])
    ylabel(sets{s})
    set(gca, "Box", "off", "XTick", 0:90:360)
    l = legend(b, sprintf('DS & OS (%d)', sum(dirTuned & oriTuned)), ...
        sprintf('DS (%d)', sum(dirTuned & ~oriTuned)));
    l.Box = "off";
    xlim([-20 380])
    xlabel('Direction (deg)')
    io.saveFigure(gcf, fPlots, sprintf('tuning_%s_directionPrefHist', sets{s}));

    % orientation
    figure
    n1 = histcounts(oriPreferences(oriTuned & dirTuned), oriEdges);
    n2 = histcounts(oriPreferences(oriTuned & ~dirTuned), oriEdges);
    % merge counts in first and last bars (around 0 and 180 deg)
    n1(1) = n1(1) + n1(end);
    n1(end) = n1(1);
    n2(1) = n2(1) + n2(end);
    n2(end) = n2(1);
    b = bar(oriEdges(1:end-1)+7.5, [n1' n2'], 'stacked');
    b(1).FaceColor = 'k';
    b(2).FaceColor = [.75 .75 .75];
    ylim([0 400])
    ylabel(sets{s})
    set(gca, "Box", "off", "XTick", 0:45:180)
    l = legend(b, sprintf('OS & DS (%d)', sum(dirTuned & oriTuned)), ...
        sprintf('OS (%d)', sum(~dirTuned & oriTuned)));
    l.Box = "off";
    xlim([-10 190])
    xlabel('Orientation (deg)')
    ylabel(sprintf('Proportion %s', sets{s}))
    io.saveFigure(gcf, fPlots, sprintf('tuning_%s_orientationPrefHist', sets{s}));

    %% Direction vs orientation preference scatterplot
    ind = dirTuned & oriTuned;
    figure
    hold on
    plot([0 180], [0 180], 'Color', [1 1 1].*0.5)
    plot([180 360], [0 180], 'Color', [1 1 1].*0.5)
    scatter(dirPreferences(ind), oriPreferences(ind), 15, 'k', 'filled')
    scatter(dirPreferences(indExamples), oriPreferences(indExamples), ...
        40, lines(length(indExamples)), 'filled')
    axis equal
    xlim([-10 370])
    ylim([-10 190])
    set(gca, "Box", "off", "XTick", 0:90:360, "YTick", 0:90:180)
    xlabel('Direction (deg)')
    ylabel('Orientation (deg)')
    title(sprintf('%s (n=%d)', sets{s}, sum(ind)))
    io.saveFigure(gcf, fPlots, sprintf('tuning_%s_dirVsOriScatter', sets{s}));

    %% DS vs OS scatterplot
    ind = dirTuned | oriTuned;
    figure
    h = gscatter(dirSel(ind), oriSel(ind), dirTuned(ind) + ...
        2*oriTuned(ind), repmat([0.4 0.8 0]', 1, 3), [], 15);
    hold on
    scatter(dirSel(indExamples), oriSel(indExamples), 40, ...
        lines(length(indExamples)), "filled")
    l = legend(h, 'DS', 'OS', 'DS & OS', "Location", "bestoutside");
    l.Box = "off";
    axis padded equal
    mini = -0.05;
    maxi = 0.85;
    axis([mini maxi mini maxi])
    set(gca, "Box", "off")
    xlabel('Direction selectivity')
    ylabel('Orientation selectivity')
    title(sprintf('%s', sets{s}))
    io.saveFigure(gcf, fPlots, ...
        sprintf('tuning_%s_selectivityDirVsOri', sets{s}));

    %% Direct comparison: boutons vs neurons + mean across datasets
    nDir = histcounts(dataset(dirTuned), 0.5:max(dataset)+1);
    nOri = histcounts(dataset(oriTuned), 0.5:max(dataset)+1);
    indSetsDir = nDir >= minUnits;
    indSetsOri = nOri >= minUnits;
    numDir(s) = sum(indSetsDir);
    numOri(s) = sum(indSetsOri);

    %% Mean preference histograms across all datasets
    % determine histograms for each dataset
    dirHists = NaN(length(dirBinsCoarse), max(dataset));
    oriHists = NaN(length(oriBins), max(dataset));
    for c = 1:max(dataset)
        pref = dirPreferences(dirTuned & c==dataset);
        dirHists(:,c) = histcounts(pref, dirEdges);
        pref = oriPreferences(oriTuned & c==dataset);
        oriHists(:,c) = histcounts(pref, oriEdges);
    end
    % merge counts in first and last bars (around 0 and 360 or 180 deg)
    dirHists(1,:) = dirHists(1,:) + dirHists(end,:);
    dirHists(end,:) = dirHists(1,:);
    oriHists(1,:) = oriHists(1,:) + oriHists(end,:);
    oriHists(end,:) = oriHists(1,:);
    % normalize histograms (to sum 1, without last bin)
    dirHists = dirHists ./ sum(dirHists(1:end-1,:),1);
    oriHists = oriHists ./ sum(oriHists(1:end-1,:),1);
    % interpolate/smooth histograms
    dirHistsSmooth = interp1([-30 dirBinsCoarse 390], ...
        dirHists([end-1 1:end 2],:), dirBinsFine, 'pchip');
    oriHistsSmooth = interp1([-15 oriBins 195], ...
        oriHists([end-1 1:end 2],:), oriBinsFine, 'pchip');

    % plot direction
    figure(figDir)
    m = mean(dirHists(:,indSetsDir),2);
    mSm = mean(dirHistsSmooth(:,indSetsDir),2);
    sSm = std(dirHistsSmooth(:,indSetsDir),0,2) ./ sqrt(sum(indSetsDir));
    fill(dirBinsFine([1:end end:-1:1]), [mSm-sSm; flip(mSm+sSm)], 'k', ...
        "FaceColor", colSets(s,:), "FaceAlpha", 0.5, "EdgeColor", "none")
    hDir(s) = plot(dirBinsFine, mSm, 'Color', colSets(s,:), "LineWidth", 1);
    plot(dirBinsCoarse, m, '.', 'Color', colSets(s,:), "MarkerSize", 30)

    % plot orientation
    figure(figOri)
    m = mean(oriHists(:,indSetsOri),2);
    mSm = mean(oriHistsSmooth(:,indSetsOri),2);
    sSm = std(oriHistsSmooth(:,indSetsOri),0,2) ./ sqrt(sum(indSetsOri));
    fill(oriBinsFine([1:end end:-1:1]), [mSm-sSm; flip(mSm+sSm)], 'k', ...
        "FaceColor", colSets(s,:), "FaceAlpha", 0.5, "EdgeColor", "none")
    hOri(s) = plot(oriBinsFine, mSm, 'Color', colSets(s,:), "LineWidth", 1);
    plot(oriBins, m, '.', 'Color', colSets(s,:), "MarkerSize", 30)

    %% Mean DS and OS histograms across datasets
    selEdges = 0:0.1:0.8;
    selBins = selEdges(1:end-1) + 0.05;
    dirHists = NaN(length(selBins), max(dataset));
    medianDS{s} = NaN(max(dataset),1);
    oriHists = NaN(length(selBins), max(dataset));
    medianOS{s} = NaN(max(dataset),1);
    for c = 1:max(dataset)
        sel = dirSel(dirTuned & c==dataset);
        dirHists(:,c) = histcounts(sel, selEdges);
        medianDS{s}(c) = median(sel);
        sel = oriSel(oriTuned & c==dataset);
        oriHists(:,c) = histcounts(sel, selEdges);
        medianOS{s}(c) = median(sel);
    end
    % normalize histograms (to sum 1)
    dirHists = dirHists ./ sum(dirHists,1);
    oriHists = oriHists ./ sum(oriHists,1);
    
    % plot direction
    figure(figDS)
    m = mean(dirHists(:,indSetsDir),2);
    sem = std(dirHists(:,indSetsDir),0,2) ./ sqrt(sum(indSetsDir));
    fill(selBins([1:end end:-1:1]), [m-sem; flip(m+sem)], 'k', ...
        "FaceColor", colSets(s,:), "FaceAlpha", 0.5, "EdgeColor", "none")
    plot(selBins, m, '.-', 'Color', colSets(s,:), "LineWidth", 1, ...
        "MarkerSize", 30)
    m = median(medianDS{s}, "omitnan");
    plot(m, 0.5, 'v', "MarkerFaceColor", colSets(s,:), ...
        "MarkerEdgeColor", "none", "MarkerSize", 10)
    
    % plot orientation
    figure(figOS)
    m = mean(oriHists(:,indSetsOri),2);
    sem = std(oriHists(:,indSetsOri),0,2) ./ sqrt(sum(indSetsOri));
    fill(selBins([1:end end:-1:1]), [m-sem; flip(m+sem)], 'k', ...
        "FaceColor", colSets(s,:), "FaceAlpha", 0.5, "EdgeColor", "none")
    plot(selBins, m, '.-', 'Color', colSets(s,:), "LineWidth", 1, ...
        "MarkerSize", 30)
    m = median(medianOS{s}, "omitnan");
    plot(m, 0.5, 'v', "MarkerFaceColor", colSets(s,:), ...
        "MarkerEdgeColor", "none", "MarkerSize", 10)

    %% Extra plots
    % plot DS and OS against brain depth
    figure('Position', [680 60 440 915])
    scatter(dirSel(dirTuned), depth(dirTuned), 'k', 'filled')
    set(gca, 'YDir', 'reverse')
    xlim([0 0.8])
    ylim([0 100])
    xlabel('DS')
    ylabel('Depth from surface')
    title(sets{s})
    io.saveFigure(gcf, fullfile(fPlots, 'extra'), ...
        sprintf('tuning_%s_DSvsDepth', sets{s}));
    figure('Position', [1150 60 440 915])
    scatter(oriSel(oriTuned), depth(oriTuned), 'k', 'filled')
    set(gca, 'YDir', 'reverse')
    xlim([0 0.8])
    ylim([0 100])
    xlabel('OS')
    ylabel('Depth from surface')
    title(sets{s})
    io.saveFigure(gcf, fullfile(fPlots, 'extra'), ...
        sprintf('tuning_%s_OSvsDepth', sets{s}));

    % plot direction/orientation preferences separately for inhibitory and
    % excitatory neurons
    if s == 2
        % plot direction preference histogram
        % 1. only inhibitory neurons
        figure
        n1 = histcounts(dirPreferences(dirTuned & oriTuned & isGad==1), dirEdges);
        n2 = histcounts(dirPreferences(dirTuned & ~oriTuned & isGad==1), dirEdges);
        b = bar(dirBinsCoarse, [n1' n2'], 'stacked');
        b(1).FaceColor = 'k';
        b(2).FaceColor = [.5 .5 .5];
        ylabel(sets{s})
        set(gca, "Box", "off", "XTick", 0:90:360)
        l = legend(b, sprintf('DS & OS (%d)', sum(n1)), ...
            sprintf('DS (%d)', sum(n2)));
        l.Box = "off";
        xlim([-10 370])
        xlabel('Direction (deg)')
        title(sprintf('Inhibitory neurons (n = %d)', sum(n1+n2)))
        io.saveFigure(gcf, fullfile(fPlots, 'extra'), ...
            sprintf('tuning_%s_inhibitory_directionPrefHist', sets{s}));
        % 2. only excitatory neurons
        figure
        n1 = histcounts(dirPreferences(dirTuned & oriTuned & isGad==-1), dirEdges);
        n2 = histcounts(dirPreferences(dirTuned & ~oriTuned & isGad==-1), dirEdges);
        b = bar(dirBinsCoarse, [n1' n2'], 'stacked');
        b(1).FaceColor = 'k';
        b(2).FaceColor = [.5 .5 .5];
        ylabel(sets{s})
        set(gca, "Box", "off", "XTick", 0:90:360)
        l = legend(b, sprintf('DS & OS (%d)', sum(n1)), ...
            sprintf('DS (%d)', sum(n2)));
        l.Box = "off";
        xlim([-10 370])
        xlabel('Direction (deg)')
        title(sprintf('Excitatory neurons (n = %d)', sum(n1+n2)))
        io.saveFigure(gcf, fullfile(fPlots, 'extra'), ...
            sprintf('tuning_%s_excitatory_directionPrefHist', sets{s}));

        % plot orientation preference histogram
        % 1. only inhibitory neurons
        figure
        n1 = histcounts(oriPreferences(dirTuned & oriTuned & isGad==1), oriEdges);
        n2 = histcounts(oriPreferences(~dirTuned & oriTuned & isGad==1), oriEdges);
        b = bar(oriBinsCoarse, [n1' n2'], 'stacked');
        b(1).FaceColor = 'k';
        b(2).FaceColor = [.8 .8 .8];
        ylabel(sets{s})
        set(gca, "Box", "off", "XTick", 0:90:180)
        l = legend(b, sprintf('DS & OS (%d)', sum(n1)), ...
            sprintf('OS (%d)', sum(n2)));
        l.Box = "off";
        xlim([-10 190])
        xlabel('Orientation (deg)')
        title(sprintf('Inhibitory neurons (n = %d)', sum(n1+n2)))
        io.saveFigure(gcf, fullfile(fPlots, 'extra'), ...
            sprintf('tuning_%s_inhibitory_orientationPrefHist', sets{s}));
        % 2. only excitatory neurons
        figure
        n1 = histcounts(oriPreferences(dirTuned & oriTuned & isGad==-1), oriEdges);
        n2 = histcounts(oriPreferences(~dirTuned & oriTuned & isGad==-1), oriEdges);
        b = bar(oriBinsCoarse, [n1' n2'], 'stacked');
        b(1).FaceColor = 'k';
        b(2).FaceColor = [.8 .8 .8];
        ylabel(sets{s})
        set(gca, "Box", "off", "XTick", 0:90:180)
        l = legend(b, sprintf('DS & OS (%d)', sum(n1)), ...
            sprintf('OS (%d)', sum(n2)));
        l.Box = "off";
        xlim([-10 190])
        xlabel('Orientation (deg)')
        title(sprintf('Excitatory neurons (n = %d)', sum(n1+n2)))
        io.saveFigure(gcf, fullfile(fPlots, 'extra'), ...
            sprintf('tuning_%s_excitatory_orientationPrefHist', sets{s}));
    end
end

figure(figDir)
set(gca, "Box", "off", "XTick", 0:90:360)
l = legend(hDir, sprintf('boutons (%d)', numDir(1)), ...
    sprintf('neurons (%d)', numDir(2)));
l.Box = "off";
xlim([-10 370])
ylim([0 0.26])
xlabel('Direction (deg)')
ylabel('Proportion of units')
io.saveFigure(gcf, fPlots, 'tuning_directionPrefHist');

figure(figOri)
set(gca, "Box", "off", "XTick", 0:45:180)
l = legend(hOri, sprintf('boutons (%d)', numOri(1)), ...
    sprintf('neurons (%d)', numOri(2)));
l.Box = "off";
xlim([-10 190])
ylim([0 0.26])
xlabel('Orientation (deg)')
ylabel('Proportion of units')
io.saveFigure(gcf, fPlots, 'tuning_orientationPrefHist');

figure(figDS)
set(gca, "Box", "off", "XTick", selEdges(1:2:end))
xlim(selEdges([1 end]))
ylim([0 0.5])
xlabel('Direction selectivity')
ylabel('Proportion of units')
io.saveFigure(gcf, fPlots, 'tuning_directionSelectivity');

figure(figOS)
set(gca, "Box", "off", "XTick", selEdges(1:2:end))
xlim(selEdges([1 end]))
ylim([0 0.5])
xlabel('Orientation selectivity')
ylabel('Proportion of units')
io.saveFigure(gcf, fPlots, 'tuning_orientationSelectivity');

%% Tests
% DS
m = [median(medianDS{1}, "omitnan"), median(medianDS{2}, "omitnan")];
fprintf('Direction selectivity (median):\n')
fprintf(' Boutons: %.4f\n', m(1))
fprintf(' Neurons: %.4f\n', m(2))
fprintf(' Comparison: boutons - neurons = %.4f (p = %.4f)\n', ...
    m(1) - m(2), ranksum(medianDS{1}, medianDS{2}))
% OS
m = [median(medianOS{1}, "omitnan"), median(medianOS{2}, "omitnan")];
fprintf('Direction selectivity (median):\n')
fprintf(' Boutons: %.4f\n', m(1))
fprintf(' Neurons: %.4f\n', m(2))
fprintf(' Comparison: boutons - neurons = %.4f (p = %.4f)\n', ...
    m(1) - m(2), ranksum(medianOS{1}, medianOS{2}))