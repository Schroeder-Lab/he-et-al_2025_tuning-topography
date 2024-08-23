%% Folders
getFolders;

%% Parameters
sets = {'boutons', 'neurons'};
measures = {'direction', 'orientation'};
retinotopyRF = true;
% retinotopyRF = false;
minEV = [0.01 0];
minPeak = [6 12];
maxP = 0.05;
gridDist = 2;
gridRadius = 5;
threshTransparent = 20;
numPerm = 1000;

%% Add paths
addpath(genpath(fullfile(folders.tools, 'npy-matlab')))
addpath(fullfile(folders.repo))

%% For all plots
fPlot = fullfile(folders.plots, 'Figure04');
if ~isfolder(fPlot)
    mkdir(fPlot)
end

%% Load data: RF position, tuning preferences
% TODO: remove units with RFs that are too close to monitor edge?

% data: .rfPos, .oriPref, .osi, .dirPref, .dsi, .set
data = struct('rfPos', cell(2,1), 'dirPref', [], 'DSI', [], ...
    'oriPref', [], 'OSI', [], 'set', []);
for s = 1:2
    subjDirs = dir(fullfile(folders.data, sets{s}, 'SS*'));
    count = 1;
    for subj = 1:length(subjDirs) % animals
        name = subjDirs(subj).name;
        dateDirs = dir(fullfile(folders.data, sets{s}, name, '2*'));
        for d = 1:length(dateDirs) %dates
            date = dateDirs(d).name;
            f = fullfile(folders.data, sets{s}, name, date);
            % ignore session if stimulus was not presented
            if ~isfile(fullfile(f, '_ss_gratingsDrifting.intervals.npy')) || ...
                    ~isfile(fullfile(f, '_ss_sparseNoise.times.npy'))
                continue
            end
            % load data
            if ~retinotopyRF
                dt = io.getRFFits(f);
                pos = dt.fitParameters(:,[2 4]);
                ind = (dt.EV > minEV(1) & dt.peaks > minPeak(1)) | ...
                    (dt.EV > minEV(2) & dt.peaks > minPeak(2));
                pos(~ind,:) = NaN;
                clear dt
            else
                if ~isfile(fullfile(f, '_ss_rf.posRetinotopy.npy'))
                    continue
                end
                pos = readNPY(fullfile(f, '_ss_rf.posRetinotopy.npy'));
                edges = readNPY(fullfile(f, '_ss_rfDescr.edges.npy'));
                ind = pos(:,1)<edges(1) | pos(:,1)>edges(2) | ...
                    pos(:,2)>edges(3) | pos(:,2)<edges(4);
                pos(ind,:) = NaN;
            end
            data(s).rfPos = [data(s).rfPos; pos];
            [dirTuning, oriTuning] = io.getTuningResults(f, 'gratingsDrifting');
            dp = dirTuning.preference;
            dsi = dirTuning.selectivity;
            dp(dirTuning.pValue >= maxP) = NaN;
            dsi(dirTuning.pValue >= maxP) = NaN;
            op = oriTuning.preference;
            osi = oriTuning.selectivity;
            op(oriTuning.pValue >= maxP) = NaN;
            osi(oriTuning.pValue >= maxP) = NaN;
            data(s).dirPref = [data(s).dirPref; dp];
            data(s).DSI = [data(s).DSI; dsi];
            data(s).oriPref = [data(s).oriPref; op];
            data(s).OSI = [data(s).OSI; osi];
            data(s).set = [data(s).set; ones(size(dp)) .* count];

            count = count + 1;
        end
    end
end

%% Scatterplot showing preferred direction/orientation of each unit
limits = [-135 -85 -17 42];
if retinotopyRF
    str = 'retinoptopic';
else
    str = 'measured';
end
for s = 1:2
    % direction
    figure
    ind = ~any(isnan(data(s).rfPos),2) & ~isnan(data(s).dirPref);
    scatter(data(s).rfPos(ind,1), data(s).rfPos(ind,2), [], ...
        data(s).dirPref(ind), "filled")
    clim([0 360])
    colormap hsv
    c = colorbar;
    c.Label.String = 'Preferred direction';
    c.Limits = [2 360];
    c.Ticks = 0:45:360;
    c.Box = "off";
    axis image
    axis(limits)
    xlabel('Azimuth (deg)')
    ylabel('Elevation (deg)')
    title(sprintf('Preferred directions of %s @ %s RFs (n=%d)', ...
        sets{s}, str, sum(ind)))
    io.saveFigure(gcf, fPlot, sprintf('globalScatter_%s_direction_%sRFs', ...
        sets{s}, str))

    % orientation
    figure
    ind = ~any(isnan(data(s).rfPos),2) & ~isnan(data(s).dirPref);
    scatter(data(s).rfPos(ind,1), data(s).rfPos(ind,2), [], ...
        data(s).oriPref(ind), "filled")
    clim([0 180])
    colormap hsv
    c = colorbar;
    c.Label.String = 'Preferred Orientation';
    c.Limits = [2 180];
    c.Ticks = 0:45:180;
    c.Box = "off";
    axis image
    axis(limits)
    xlabel('Azimuth (deg)')
    ylabel('Elevation (deg)')
    title(sprintf('Preferred orientations of %s @ %s RFs (n=%d)', ...
        sets{s}, str, sum(ind)))

    io.saveFigure(gcf, fPlot, sprintf('globalScatter_%s_orientation_%sRFs', ...
        sets{s}, str))
end

%% Create smoothed maps and surrogate maps
maps = struct('x', cell(2,1), 'y', [], 'direction', [], 'orientation', []);
limits = [-134 -94 -16 40];
for s = 1:2
    % direction
    [x, y, preferences, consistencies, counts] = ...
        spatial.makeSmoothMap(data(s).rfPos, data(s).dirPref, ...
        gridDist, gridRadius, limits);
    maps(s).x = x;
    maps(s).y = y;
    maps(s).direction.preferences = preferences;
    maps(s).direction.consistencies = consistencies;
    maps(s).direction.counts = counts;
    maps(s).direction.nullPrefs = NaN([size(preferences) numPerm]);
    maps(s).direction.nullCons = NaN([size(preferences) numPerm]);

    valid = find(~any(isnan([data(s).rfPos, data(s).dirPref]), 2));
    rng('default');
    for p = 1:numPerm
        order = randperm(length(valid));
        [~, ~, prefs, cons] = spatial.makeSmoothMap(...
            data(s).rfPos(valid,:), data(s).dirPref(valid(order)), ...
            gridDist, gridRadius, limits);
        maps(s).direction.nullPrefs(:,:,p) = prefs;
        maps(s).direction.nullCons(:,:,p) = cons;
    end

    % orientation
    [x, y, preferences, consistencies, counts] = ...
        spatial.makeSmoothMap(data(s).rfPos, data(s).oriPref, ...
        gridDist, gridRadius, limits);
    maps(s).orientation.preferences = preferences;
    maps(s).orientation.consistencies = consistencies;
    maps(s).orientation.counts = counts;
    maps(s).orientation.nullPrefs = NaN([size(preferences) numPerm]);
    maps(s).orientation.nullCons = NaN([size(preferences) numPerm]);

    valid = find(~any(isnan([data(s).rfPos, data(s).oriPref]), 2));
    rng('default');
    for p = 1:numPerm
        order = randperm(length(valid));
        [~, ~, prefs, cons] = spatial.makeSmoothMap(...
            data(s).rfPos(valid,:), data(s).oriPref(valid(order)), ...
            gridDist, gridRadius, limits);
        maps(s).orientation.nullPrefs(:,:,p) = prefs;
        maps(s).orientation.nullCons(:,:,p) = cons;
    end
end

%% Plot smoothed preference maps pooling datasets
if retinotopyRF
    str = 'retinoptopic';
else
    str = 'measured';
end
for s = 1:2
    % direction
    transparency = maps(s).direction.counts./threshTransparent;
    transparency(transparency>1) = 1;
    % preferences
    figure
    imagesc(x([1 end]), y([1 end]), maps(s).direction.preferences, ...
        'AlphaData', transparency, [0 360])
    colormap([[1 1 1]; hsv])
    c = colorbar;
    c.Label.String = 'Preferred Direction';
    c.Limits = [2 360];
    c.Ticks = 0:45:360;
    c.Box = "off";
    axis image
    set(gca, "Box", "off", "YDir", "normal")
    xlabel('Azimuth (deg)')
    ylabel('Elevation (deg)')
    title(sprintf('Global direction map for %s (%s RFs)', sets{s}, str))
    io.saveFigure(gcf, fPlot, sprintf('globalMap_%s_direction_%sRFs', ...
        sets{s}, str))
    % consistency
    figure
    cons = maps(s).direction.consistencies;
    cons(isnan(cons)) = 0;
    imagesc(x([1 end]), y([1 end]), cons)
    colormap(gray)
    c = colorbar;
    c.Label.String = 'Consistency';
    c.Limits = [0 1];
    c.Box = "off";
    axis image
    set(gca, "Box", "off", "YDir", "normal")
    xlabel('Azimuth (deg)')
    ylabel('Elevation (deg)')
    title(sprintf('Global direction consistency map for %s (%s RFs)', sets{s}, str))
    io.saveFigure(gcf, fPlot, sprintf('globalMap_%s_consistency-of-direction_%sRFs', ...
        sets{s}, str))
    % preference * consistency
    colors = hsv(360);
    prefs = round(maps(s).direction.preferences(:));
    prefs(prefs==0) = 360;
    ind = ~isnan(prefs);
    im = ones([size(prefs), 3]) .* 0.5;
    im(ind,1,:) = colors(prefs(ind),:);
    im(ind,1,:) = maps(s).direction.consistencies(ind) .* im(ind,1,:) + ...
        (1-maps(s).direction.consistencies(ind)) .* ...
        ones(size(prefs(ind))) .* 0.5;
    im = reshape(im, size(maps(s).direction.preferences,1), [], 3);
    figure
    image(x([1 end]), y([1 end]), im)
    axis image
    set(gca, "Box", "off", "YDir", "normal")
    xlabel('Azimuth (deg)')
    ylabel('Elevation (deg)')
    title(sprintf('Global direction * consistency map for %s (%s RFs)', sets{s}, str))
    io.saveFigure(gcf, fPlot, sprintf('globalMap_%s_direction-x-consistency_%sRFs', ...
        sets{s}, str))

    % orientation
    transparency = maps(s).orientation.counts./threshTransparent;
    transparency(transparency>1) = 1;
    % preferences
    figure
    imagesc(x([1 end]), y([1 end]), maps(s).orientation.preferences,...
        'AlphaData', transparency, [0 180])
    colormap([[1 1 1]; hsv])
    c = colorbar;
    c.Label.String = 'Preferred Orientation';
    c.Limits = [2 180];
    c.Ticks = 0:45:180;
    c.Box = "off";
    axis image
    set(gca, "Box", "off", "YDir", "normal")
    xlabel('Azimuth (deg)')
    ylabel('Elevation (deg)')
    title(sprintf('Global orientation map for %s (%s RFs)', sets{s}, str))
    io.saveFigure(gcf, fPlot, sprintf('globalMap_%s_orientation_%sRFs', ...
        sets{s}, str))
    % consistency
    figure
    cons = maps(s).orientation.consistencies;
    cons(isnan(cons)) = 0;
    imagesc(x([1 end]), y([1 end]), cons)
    colormap(gray)
    c = colorbar;
    c.Label.String = 'Consistency';
    c.Limits = [0 1];
    c.Box = "off";
    axis image
    set(gca, "Box", "off", "YDir", "normal")
    xlabel('Azimuth (deg)')
    ylabel('Elevation (deg)')
    title(sprintf('Global orientation consistency map for %s (%s RFs)', sets{s}, str))
    io.saveFigure(gcf, fPlot, sprintf('globalMap_%s_consistency-of-orientation_%sRFs', ...
        sets{s}, str))
    % preference * consistency
    colors = hsv(180);
    prefs = round(maps(s).orientation.preferences(:));
    prefs(prefs==0) = 180;
    ind = ~isnan(prefs);
    im = ones([size(prefs), 3]) .* 0.5;
    im(ind,1,:) = colors(round(prefs(ind)),:);
    im(ind,1,:) = maps(s).orientation.consistencies(ind) .* im(ind,1,:) + ...
        (1-maps(s).orientation.consistencies(ind)) .* ...
        ones(size(prefs(ind))) .* 0.5;
    im = reshape(im, size(maps(s).orientation.preferences,1), [], 3);
    figure
    image(x([1 end]), y([1 end]), im)
    axis image
    set(gca, "Box", "off", "YDir", "normal")
    xlabel('Azimuth (deg)')
    ylabel('Elevation (deg)')
    title(sprintf('Global orientation * consistency map for %s (%s RFs)', sets{s}, str))
    io.saveFigure(gcf, fPlot, sprintf('globalMap_%s_orientation-x-consistency_%sRFs', ...
        sets{s}, str))
end

%% Plot histograms + scatters: consistencies compared to null distribution
binEdges = 0:0.05:1;
bins = binEdges(2:end) - 0.025;
binsSmooth = linspace(bins(1), bins(end), 100);
bins2 = 0.01:0.02:1;
if retinotopyRF
    str = 'retinoptopic';
else
    str = 'measured';
end
for s = 1:2
    for m = 1:2
        n = histcounts(maps(s).(measures{m}).consistencies(:), binEdges);
        n = interp1(bins, n, binsSmooth, "pchip");
        n = n ./ sum(n);
        nulls = NaN(length(bins), numPerm);
        for p = 1:numPerm
            nulls(:,p) = histcounts(maps(s).(measures{m}).nullCons(:,:,p), ...
                binEdges);
        end
        confIntv = prctile(nulls, [2.5 50 97.5], 2);
        confIntv = interp1(bins, confIntv, binsSmooth, "pchip");
        confIntv = confIntv ./ sum(confIntv(:,2));

        % histogram
        h = [0 0];
        figure
        hold on
        fill([binsSmooth flip(binsSmooth)], ...
            [confIntv(:,1); flip(confIntv(:,3))], ...
            'k', "FaceColor", 'k', "FaceAlpha", 0.3, "EdgeColor", "none")
        h(1) = plot(binsSmooth, confIntv(:,2), 'k', "LineWidth", 1);
        h(2) = plot(binsSmooth, n, 'r', "LineWidth", 1);
        legend(h, 'null', 'original')
        xlabel('Consistency')
        ylabel('Probability')
        title(sprintf('%s consistency (%s RFs) - %s', ...
            measures{m}, str, sets{s}))
        io.saveFigure(gcf, fPlot, sprintf('consistency_%s_%s_histogram_%sRFs', ...
            sets{s}, measures{m}, str))

        % scatter: consistency vs unit count
        [x, y] = meshgrid(bins2, ...
            2.5:5:max(maps(s).(measures{m}).counts,[],"all")+5);
        f = ksdensity([maps(s).(measures{m}).nullCons(:), ...
            repmat(maps(s).(measures{m}).counts(:),numPerm,1)], [x(:) y(:)]);
        f = reshape(f, size(x));
        figure
        imagesc(bins2([1 end]), y([1 end]), f)
        colormap(flipud(gray))
        set(gca, "Box", "off", "YDir", "normal")
        hold on
        scatter(maps(s).(measures{m}).consistencies, ...
            maps(s).(measures{m}).counts, 15, 'r', 'filled')
        xlabel('Consistency')
        ylabel('#units per patch')
        title(sprintf('%s consistency (%s RFs) - %s', measures{m}, str, sets{s}))
        io.saveFigure(gcf, fPlot, sprintf('consistency_%s_%s_scatter-count_%sRFs', ...
            sets{s}, measures{m}, str))
    end
end

%% Determine consistencies and null distributions for each dataset
consistency = struct('direction', cell(2,1), 'orientation', []);
strPrefs = {'dirPref', 'oriPref'};
for s = 1:2
    numRecs = max(data(s).set);
    consistency(s).direction.original = cell(numRecs,1);
    consistency(s).direction.null = cell(numRecs,1);
    consistency(s).orientation.original = cell(numRecs,1);
    consistency(s).orientation.null = cell(numRecs,1);
    for d = 1:numRecs
        indRec = data(s).set == d;
        for m = 1:2
            prefs = data(s).(strPrefs{m})(indRec);
            rfs = data(s).rfPos(indRec,:);
            [x,y,~,consistencies] = ...
                spatial.makeSmoothMap(rfs, prefs, gridDist, gridRadius);
            % only consider datasets where RF positions span a distance of
            % at least gridRadius (otherwhise the permuted null data will
            % give the same results as the original data)
            if isempty(x) || diff(x([1 end]))<gridRadius || ...
                    diff(y([1 end]))<gridRadius
                continue
            end
            consistency(s).(measures{m}).original{d} = consistencies;
            consistency(s).(measures{m}).null{d} = ...
                NaN([size(consistencies) numPerm]);

            valid = find(~any(isnan([rfs, prefs]), 2));
            rng('default');
            for p = 1:numPerm
                order = randperm(length(valid));
                [~, ~,~, cons] = spatial.makeSmoothMap(...
                    rfs(valid,:), prefs(valid(order)), ...
                    gridDist, gridRadius);
                consistency(s).(measures{m}).null{d}(:,:,p) = cons;
            end
        end
    end
end

%% Plot consistencies compared to null distribution, per dataset
binEdges = 0:0.1:1;
bins = binEdges(2:end) - 0.05;
binsSmooth = linspace(bins(1), bins(end), 100);
if retinotopyRF
    str = 'retinoptopic';
else
    str = 'measured';
end
for s = 1:2
    for m = 1:2
        probs = NaN(length(bins), ...
            length(consistency(s).(measures{m}).original));
        probsNull = NaN(length(bins), ...
            length(consistency(s).(measures{m}).original));
        probsNullDistr = NaN(length(bins), ...
            length(consistency(s).(measures{m}).original), numPerm);
        for d = 1:length(consistency(s).(measures{m}).original)
            if isempty(consistency(s).(measures{m}).original{d})
                continue
            end
            probs(:,d) = histcounts( ...
                consistency(s).(measures{m}).original{d}, binEdges, ...
                "Normalization", "probability");
            probsNull(:,d) = histcounts( ...
                consistency(s).(measures{m}).null{d}, binEdges, ...
                "Normalization", "probability");
            for p = 1:numPerm
                probsNullDistr(:,d,p) = histcounts( ...
                consistency(s).(measures{m}).null{d}(:,:,p), binEdges, ...
                "Normalization", "probability");
            end
        end
        valid = ~any(isnan(probs),1);
        lineColors = turbo(length(valid));

        % consistencies of null data
        smNull = NaN(length(binsSmooth), length(valid));
        smNull(:,valid) = interp1(bins, probsNull(:,valid), binsSmooth, "pchip");
        smNull = smNull ./ sum(smNull,1);
        % consistencies of original data
        smOrig = NaN(length(binsSmooth), length(valid));
        smOrig(:,valid) = interp1(bins, probs(:,valid), binsSmooth, "pchip");
        smOrig = smOrig ./ sum(smOrig,1);
        figure
        hold on
        for v = 1:length(valid)
            if ~valid(v)
                continue
            end
            % consistencies of null data
            plot(binsSmooth, smNull(:,v), ':', "LineWidth", 1, ...
                "Color", lineColors(v,:));
            % consistencies of original data
            plot(binsSmooth, smOrig(:,v), "LineWidth", 1, ...
                "Color", lineColors(v,:));
        end
        h = [0 0];
        h(1) = plot(binsSmooth, mean(smNull,2,"omitnan"), 'k:', ...
            "LineWidth", 2);
        h(2) = plot(binsSmooth, mean(smOrig,2,"omitnan"), 'k', ...
            "LineWidth", 2);
        set(gca, "Box", "off")
        legend(h, 'null', 'original')
        xlabel(sprintf('Consistency of %s preferences', measures{m}))
        ylabel('Probability')
        title(sprintf('%s consistency (%sRFs) - %s', ...
            measures{m}, str, sets{s}))
        io.saveFigure(gcf, fPlot, sprintf(...
            'consistency_%s_%s_histograms-per-dataset_%sRFs', ...
            sets{s}, measures{m}, str))

        % mean consistencies compared to null distribution of mean
        mn = cellfun(@mean, consistency(s).(measures{m}).original, ...
            repmat({"all"},length(valid),1), ...
            repmat({"omitnan"},length(valid),1));
        nullM = cellfun(@mean, consistency(s).(measures{m}).null, ...
            repmat({[1 2]},length(valid),1), ...
            repmat({"omitnan"},length(valid),1), 'UniformOutput', false);
        nullM = cellfun(@squeeze, nullM, 'UniformOutput', false);
        nullM_flat = NaN(numPerm, length(valid));
        nullM_flat(:,valid) = cat(2, nullM{valid});
        figure
        hold on
        h = [0 0 0];
        for v = 1:length(valid)
            if ~valid(v)
                continue
            end
            inv = prctile(nullM_flat(:,v), [2.5 97.5]);
            h(1) = plot(inv, [1 1].*v, ...
                "Color", lineColors(v,:));
            h(2) = plot(prctile(nullM_flat(:,v), 50), v, 'o', ...
                "MarkerFaceColor", lineColors(v,:), "MarkerEdgeColor", "none");
            h(3) = plot(mn(v), v, 'v', ...
                "MarkerFaceColor", lineColors(v,:), "MarkerEdgeColor", "none");
            if mn(v)>inv(2) || mn(v)<inv(1)
                plot(1.1, v, 'k*')
            end
        end
        set(gca, "Box", "off")
        legend(h, 'null interval', 'null median', 'original', "Location", "bestoutside")
        xlim([0 1.15])
        ylim([find(valid,1)-1 find(valid,1,'last')+1])
        xlabel(sprintf('Mean consistency of %s preferences', measures{m}))
        ylabel('Dataset')
        title('Original relative to null')
        io.saveFigure(gcf, fPlot, sprintf(...
            'consistency_%s_%s_intervals-per-dataset_%sRFs', ...
            sets{s}, measures{m}, str))
    end
end