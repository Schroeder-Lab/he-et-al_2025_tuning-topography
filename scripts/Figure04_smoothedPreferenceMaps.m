function maps = Figure05_smoothedPreferenceMaps(glob, fPlots, data, ...
    sets, retinotopyRF)

%% Parameters
% create surrogate global maps of direction/orientation tuning
numPerm = 1000;

% plotting global maps of direction/orientation tuning
gridDist = 2;
patchRadius = 5;
threshTransparent = 20;
colorsDir = colmaps.colorcet('C7'); % has 256 entries/colours
colorsOri = colmaps.colorcet('C1'); % has 256 entries/colours

%% Create smoothed maps and surrogate maps
maps = struct('x', cell(2,1), 'y', [], 'direction', [], 'orientation', []);
limits = [-132 -84 -16 38];
for s = 1:2
    % direction
    [x, y, preferences, consistencies, counts] = ...
        spatial.makeSmoothMap(data(s).rfPos, data(s).dirPref, ...
        gridDist, patchRadius, limits);
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
            gridDist, patchRadius, limits);
        maps(s).direction.nullPrefs(:,:,p) = prefs;
        maps(s).direction.nullCons(:,:,p) = cons;
    end

    % orientation
    [~, ~, preferences, consistencies, counts] = ...
        spatial.makeSmoothMap(data(s).rfPos, data(s).oriPref .* 2, ...
        gridDist, patchRadius, limits);
    maps(s).orientation.preferences = preferences ./ 2;
    maps(s).orientation.consistencies = consistencies;
    maps(s).orientation.counts = counts;
    maps(s).orientation.nullPrefs = NaN([size(preferences) numPerm]);
    maps(s).orientation.nullCons = NaN([size(preferences) numPerm]);

    valid = find(~any(isnan([data(s).rfPos, data(s).oriPref]), 2));
    rng('default');
    for p = 1:numPerm
        order = randperm(length(valid));
        [~, ~, prefs, cons] = spatial.makeSmoothMap(...
            data(s).rfPos(valid,:), data(s).oriPref(valid(order)) .* 2, ...
            gridDist, patchRadius, limits);
        maps(s).orientation.nullPrefs(:,:,p) = prefs ./ 2;
        maps(s).orientation.nullCons(:,:,p) = cons;
    end
end

%% Plot smoothed preference maps pooling datasets
for s = 1:2
    if retinotopyRF(s)
        str = 'retinoptopic';
    else
        str = 'measured';
    end
    % direction
    transparency = maps(s).direction.counts./threshTransparent;
    transparency(transparency>1) = 1;
    angles = maps(s).direction.preferences;
    angles(angles < 1) = 360;
    angles = round(angles / 360 * 256);
    % preferences
    figure('Position', glob.figPositionDefault)
    imagesc(maps(s).x([1 end]), maps(s).y([1 end]), angles, ...
        'AlphaData', transparency, [0 256])
    colormap([[1 1 1]; colorsDir])
    c = colorbar;
    c.Label.String = 'Preferred Direction';
    c.Limits = [0 256];
    c.Ticks = (0:45:360) ./ 360 .* 256;
    c.TickLabels = 0:45:360;
    c.Box = "off";
    axis image
    set(gca, "Box", "off", "YDir", "normal")
    xlabel('Azimuth (deg)')
    ylabel('Elevation (deg)')
    title(sprintf('Global direction map for %s (%s RFs)', sets{s}, str))
    if ~isempty(fPlots)
        io.saveFigure(gcf, fPlots, sprintf('globalMap_%s_direction_%sRFs', ...
            sets{s}, str))
    end
    % consistency
    figure('Position', glob.figPositionDefault)
    cons = maps(s).direction.consistencies;
    cons(isnan(cons)) = 0;
    imagesc(maps(s).x([1 end]), maps(s).y([1 end]), cons, [0 1])
    colormap(flip(gray))
    c = colorbar;
    c.Label.String = 'Consistency';
    c.Limits = [0 1];
    c.Box = "off";
    axis image
    set(gca, "Box", "off", "YDir", "normal")
    xlabel('Azimuth (deg)')
    ylabel('Elevation (deg)')
    title(sprintf('Global direction consistency map for %s (%s RFs)', ...
        sets{s}, str))
    if ~isempty(fPlots)
        io.saveFigure(gcf, fPlots, sprintf('globalMap_%s_consistency-of-direction_%sRFs', ...
            sets{s}, str))
    end

    % orientation
    transparency = maps(s).orientation.counts./threshTransparent;
    transparency(transparency>1) = 1;
    angles = maps(s).orientation.preferences;
    angles(angles < 1) = 180;
    % preferences
    figure('Position', glob.figPositionDefault)
    imagesc(maps(s).x([1 end]), maps(s).y([1 end]), angles, ...
        'AlphaData', transparency, [0 180])
    colormap([[1 1 1]; colorsOri])
    c = colorbar;
    c.Label.String = 'Preferred Orientation';
    c.Limits = [0 180];
    c.Ticks = 0:45:180;
    c.Box = "off";
    axis image
    set(gca, "Box", "off", "YDir", "normal")
    xlabel('Azimuth (deg)')
    ylabel('Elevation (deg)')
    title(sprintf('Global orientation map for %s (%s RFs)', sets{s}, str))
    if ~isempty(fPlots)
        io.saveFigure(gcf, fPlots, sprintf('globalMap_%s_orientation_%sRFs', ...
            sets{s}, str))
    end
    % consistency
    figure('Position', glob.figPositionDefault)
    cons = maps(s).orientation.consistencies;
    cons(isnan(cons)) = 0;
    imagesc(maps(s).x([1 end]), maps(s).y([1 end]), cons, [0 1])
    colormap(flip(gray))
    c = colorbar;
    c.Label.String = 'Consistency';
    c.Limits = [0 1];
    c.Box = "off";
    axis image
    set(gca, "Box", "off", "YDir", "normal")
    xlabel('Azimuth (deg)')
    ylabel('Elevation (deg)')
    title(sprintf('Global orientation consistency map for %s (%s RFs)', ...
        sets{s}, str))
    if ~isempty(fPlots)
        io.saveFigure(gcf, fPlots, sprintf('globalMap_%s_consistency-of-orientation_%sRFs', ...
            sets{s}, str))
    end
end