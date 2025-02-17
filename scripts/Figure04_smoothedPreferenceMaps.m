function maps = Figure04_smoothedPreferenceMaps(data, fPlots, sets, ...
    retinotopyRF)

%% Parameters
% create surrogate global maps of direction/orientation tuning
numPerm = 1000;

% plotting global maps of direction/orientation tuning
gridDist = 2;
gridRadius = 5;
threshTransparent = 20;
colors = colmaps.colorcet('C7');

%% Create smoothed maps and surrogate maps
maps = struct('x', cell(2,1), 'y', [], 'direction', [], 'orientation', []);
limits = [-132 -84 -16 38];
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
    [~, ~, preferences, consistencies, counts] = ...
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
for s = 1:2
    if retinotopyRF(s)
        str = 'retinoptopic';
    else
        str = 'measured';
    end
    % direction
    transparency = maps(s).direction.counts./threshTransparent;
    transparency(transparency>1) = 1;
    % preferences
    figure
    imagesc(maps(s).x([1 end]), maps(s).y([1 end]), ...
        maps(s).direction.preferences, ...
        'AlphaData', transparency, [0 360])
    colormap([[1 1 1]; colors])
    c = colorbar;
    c.Label.String = 'Preferred Direction';
    c.Limits = [1 360];
    c.Ticks = [1 45:45:360];
    c.Box = "off";
    axis image
    set(gca, "Box", "off", "YDir", "normal")
    xlabel('Azimuth (deg)')
    ylabel('Elevation (deg)')
    title(sprintf('Global direction map for %s (%s RFs)', sets{s}, str))
    io.saveFigure(gcf, fPlots, sprintf('globalMap_%s_direction_%sRFs', ...
        sets{s}, str))
    % consistency
    figure
    cons = maps(s).direction.consistencies;
    cons(isnan(cons)) = 0;
    imagesc(maps(s).x([1 end]), maps(s).y([1 end]), cons)
    colormap(gray)
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
    io.saveFigure(gcf, fPlots, sprintf('globalMap_%s_consistency-of-direction_%sRFs', ...
        sets{s}, str))
    % preference * consistency
    prefs = round(maps(s).direction.preferences(:));
    prefs(prefs==0) = 360;
    % transform 360 angles to 256 colours
    prefs = round(prefs/360*256);
    ind = ~isnan(prefs);
    im = ones([size(prefs), 3]) .* 0.5;
    im(ind,1,:) = colors(prefs(ind),:);
    im(ind,1,:) = maps(s).direction.consistencies(ind) .* im(ind,1,:) + ...
        (1-maps(s).direction.consistencies(ind)) .* ...
        ones(size(prefs(ind))) .* 0.5;
    im = reshape(im, size(maps(s).direction.preferences,1), [], 3);
    figure
    image(maps(s).x([1 end]), maps(s).y([1 end]), im)
    axis image
    set(gca, "Box", "off", "YDir", "normal")
    xlabel('Azimuth (deg)')
    ylabel('Elevation (deg)')
    title(sprintf('Global direction * consistency map for %s (%s RFs)', ...
        sets{s}, str))
    io.saveFigure(gcf, fPlots, sprintf('globalMap_%s_direction-x-consistency_%sRFs', ...
        sets{s}, str))

    % orientation
    transparency = maps(s).orientation.counts./threshTransparent;
    transparency(transparency>1) = 1;
    % preferences
    figure
    imagesc(maps(s).x([1 end]), maps(s).y([1 end]), ...
        maps(s).orientation.preferences,...
        'AlphaData', transparency, [0 180])
    colormap([[1 1 1]; colors])
    c = colorbar;
    c.Label.String = 'Preferred Orientation';
    c.Limits = [1 180];
    c.Ticks = [1 45:45:180];
    c.Box = "off";
    axis image
    set(gca, "Box", "off", "YDir", "normal")
    xlabel('Azimuth (deg)')
    ylabel('Elevation (deg)')
    title(sprintf('Global orientation map for %s (%s RFs)', sets{s}, str))
    io.saveFigure(gcf, fPlots, sprintf('globalMap_%s_orientation_%sRFs', ...
        sets{s}, str))
    % consistency
    figure
    cons = maps(s).orientation.consistencies;
    cons(isnan(cons)) = 0;
    imagesc(maps(s).x([1 end]), maps(s).y([1 end]), cons)
    colormap(gray)
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
    io.saveFigure(gcf, fPlots, sprintf('globalMap_%s_consistency-of-orientation_%sRFs', ...
        sets{s}, str))
    % preference * consistency
    prefs = round(maps(s).orientation.preferences(:));
    % transform 180 angles to 256 colours
    prefs = round(prefs/180*256);
    prefs(prefs==0) = 256;
    ind = ~isnan(prefs);
    im = ones([size(prefs), 3]) .* 0.5;
    im(ind,1,:) = colors(round(prefs(ind)),:);
    im(ind,1,:) = maps(s).orientation.consistencies(ind) .* im(ind,1,:) + ...
        (1-maps(s).orientation.consistencies(ind)) .* ...
        ones(size(prefs(ind))) .* 0.5;
    im = reshape(im, size(maps(s).orientation.preferences,1), [], 3);
    figure
    image(maps(s).x([1 end]), maps(s).y([1 end]), im)
    axis image
    set(gca, "Box", "off", "YDir", "normal")
    xlabel('Azimuth (deg)')
    ylabel('Elevation (deg)')
    title(sprintf('Global orientation * consistency map for %s (%s RFs)', ...
        sets{s}, str))
    io.saveFigure(gcf, fPlots, sprintf('globalMap_%s_orientation-x-consistency_%sRFs', ...
        sets{s}, str))
end