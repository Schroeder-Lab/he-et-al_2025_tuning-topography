%% Folders
getFolders;

%% Parameters
sets = {'boutons', 'neurons'};
measures = {'direction', 'orientation'};
% retinotopyRF = true;
retinotopyRF = false;
minEV = [0.01 0];
minPeak = [6 12];
maxP = 0.05;
% TODO: introduce different grid parameters for boutons and neurons
gridDist = 2;
gridRadius = 10;
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
    axis([min(data(s).rfPos(:,1))-2 max(data(s).rfPos(:,1))+2 ...
        min(data(s).rfPos(:,2))-2 max(data(s).rfPos(:,2))+2])
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
    axis([min(data(s).rfPos(:,1))-2 max(data(s).rfPos(:,1))+2 ...
        min(data(s).rfPos(:,2))-2 max(data(s).rfPos(:,2))+2])
    xlabel('Azimuth (deg)')
    ylabel('Elevation (deg)')
    title(sprintf('Preferred orientations of %s @ %s RFs (n=%d)', ...
        sets{s}, str, sum(ind)))

    io.saveFigure(gcf, fPlot, sprintf('globalScatter_%s_orientation_%sRFs', ...
        sets{s}, str))
end

%% Create smoothed maps and surrogate maps
maps = struct('x', cell(2,1), 'y', [], 'direction', [], 'orientation', []);
for s = 1:2

    % direction
    [x, y, preferences, consistencies, counts] = ...
        spatial.makeSmoothMap(data(s).rfPos, data(s).dirPref, ...
        gridDist, gridRadius);
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
            gridDist, gridRadius);
        maps(s).direction.nullPrefs(:,:,p) = prefs;
        maps(s).direction.nullCons(:,:,p) = cons;
    end

    % orientation
    [x, y, preferences, consistencies, counts] = ...
        spatial.makeSmoothMap(data(s).rfPos, data(s).oriPref, ...
        gridDist, gridRadius);
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
            gridDist, gridRadius);
        maps(s).orientation.nullPrefs(:,:,p) = prefs;
        maps(s).orientation.nullCons(:,:,p) = cons;
    end
end

%% Smoothed preference maps pooling datasets
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

%% Histograms + scatters: consistencies compared to null distribution
binEdges = 0:0.02:1;
bins = binEdges(2:end) - 0.01;
if retinotopyRF
    str = 'retinoptopic';
else
    str = 'measured';
end
for s = 1:2
    for m = 1:2
        n = histcounts(maps(s).(measures{m}).consistencies(:), binEdges);
        nulls = NaN(length(binEdges)-1, numPerm);
        for p = 1:numPerm
            nulls(:,p) = histcounts(maps(s).(measures{m}).nullCons(:,:,p), binEdges);
        end
        confIntv = prctile(nulls, [2.5 50 97.5], 2);

        % histogram
        h = [0 0];
        figure
        hold on
        fill([bins flip(bins)], [confIntv(:,1); flip(confIntv(:,3))], ...
            'k', "FaceColor", 'k', "FaceAlpha", 0.3, "EdgeColor", "none")
        h(1) = plot(bins, confIntv(:,2), 'k', "LineWidth", 1);
        h(2) = plot(bins, n, 'r', "LineWidth", 1);
        legend(h, 'null distribution', 'real')
        xlabel('Consistency')
        ylabel('Count')
        title(sprintf('%s consistency (%s RFs) - %s', measures{m}, str, sets{s}))
        io.saveFigure(gcf, fPlot, sprintf('consistency_%s_%s_histogram_%sRFs', ...
            sets{s}, measures{m}, str))

        % scatter: consistency vs unit count
        [x, y] = meshgrid(bins, 2.5:5:max(maps(s).(measures{m}).counts,[],"all"));
        f = ksdensity([maps(s).(measures{m}).nullCons(:), ...
            repmat(maps(s).(measures{m}).counts(:),numPerm,1)], [x(:) y(:)]);
        f = reshape(f, size(x));
        figure
        imagesc(bins([1 end]), y([1 end]), f)
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
for s = 1:2
    for m = 1:2

        % CONTINUE HERE

        numRecs = max(data(s))
        for d = 1:
            % direction
            [x, y, preferences, consistencies, counts] = ...
                spatial.makeSmoothMap(data(s).rfPos, data(s).dirPref, ...
                gridDist, gridRadius);
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
                    gridDist, gridRadius);
                maps(s).direction.nullPrefs(:,:,p) = prefs;
                maps(s).direction.nullCons(:,:,p) = cons;
            end

            % orientation
            [x, y, preferences, consistencies, counts] = ...
                spatial.makeSmoothMap(data(s).rfPos, data(s).oriPref, ...
                gridDist, gridRadius);
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
                    gridDist, gridRadius);
                maps(s).orientation.nullPrefs(:,:,p) = prefs;
                maps(s).orientation.nullCons(:,:,p) = cons;
            end
        end
    end
end

%% Consistencies compared to null distribution, per dataset
