function Figure04_preferenceMapsAcrossAllDatasets(data, fPlots, sets, ...
    retinotopyRF, selectivityThresholds)

if nargin < 5 
    selectivityThresholds = [0 1; 0 1];
end

minUnits = 3;

limits = [-132 -84 -16 38];
for s = 1:2
    if retinotopyRF(s)
        str = 'retinoptopic';
    else
        str = 'measured';
    end
    % direction
    colors = colmaps.colorcet('C7'); % has 256 entries/colours
    % only use units with minimum DSI and maximum OSI
    valid = data(s).DSI >= selectivityThresholds(1,1) & ...
        data(s).OSI <= selectivityThresholds(1,2);
    % set directions to range 1-360
    angles = round(mod(data(s).dirPref,360));
    angles(angles == 0) = 360;
    figure
    ind = valid & ~any(isnan(data(s).rfPos),2) & ~isnan(data(s).dirPref);
    scatter(data(s).rfPos(ind,1), data(s).rfPos(ind,2), [], ...
        angles(ind), "filled")
    hold on
    % plot Gaussian fit (bivariate normal distribution)
    for k = 1:max(data(s).set)
        ind = data(s).set == k & ~any(isnan(data(s).rfPos),2);
        if sum(ind) < minUnits
            continue
        end
        [x, y] = algebra.getGaussianContour(data(s).rfPos(ind,1), ...
            data(s).rfPos(ind,2));
        plot(x, y, 'k')
    end
    clim([0 360])
    colormap(colors)
    c = colorbar;
    c.Label.String = 'Preferred direction';
    c.Limits = [0 360];
    c.Ticks = 0:45:360;
    c.Box = "off";
    axis image
    axis(limits)
    xlabel('Azimuth (deg)')
    ylabel('Elevation (deg)')
    n = sum(~any(isnan(data(s).rfPos),2) & ~isnan(data(s).dirPref) & valid);
    title(sprintf('%s @ %s RFs (n = %d)', sets{s}, str, n))
    io.saveFigure(gcf, fPlots, sprintf('globalScatter_%s_direction_%sRFs', ...
        sets{s}, str))

    % orientation
    colors = colmaps.colorcet('C1');
    % only use units with minimum OSI and maximum DSI
    valid = data(s).OSI >= selectivityThresholds(2,1) & ...
        data(s).DSI <= selectivityThresholds(2,2);
    % set orientations to range 1-180
    angles = round(mod(data(s).oriPref,180));
    angles(angles == 0) = 180;
    figure
    ind = valid & ~any(isnan(data(s).rfPos),2) & ~isnan(data(s).oriPref);
    scatter(data(s).rfPos(ind,1), data(s).rfPos(ind,2), [], ...
        angles(ind), "filled")
    hold on
    % plot Gaussian fit (bivariate normal distribution)
    for k = 1:max(data(s).set)
        ind = data(s).set == k & ~any(isnan(data(s).rfPos),2);
        if sum(ind) < minUnits
            continue
        end
        [x, y] = algebra.getGaussianContour(data(s).rfPos(ind,1), ...
            data(s).rfPos(ind,2));
        plot(x, y, 'k')
    end
    clim([0 180])
    colormap(colors)
    c = colorbar;
    c.Label.String = 'Preferred Orientation';
    c.Limits = [0 180];
    c.Ticks = 0:45:180;
    c.Box = "off";
    axis image
    axis(limits)
    xlabel('Azimuth (deg)')
    ylabel('Elevation (deg)')
    n = sum(~any(isnan(data(s).rfPos),2) & ~isnan(data(s).oriPref) & valid);
    title(sprintf('%s @ %s RFs (n = %d)', sets{s}, str, n))
    io.saveFigure(gcf, fPlots, sprintf('globalScatter_%s_orientation_%sRFs', ...
        sets{s}, str))
end