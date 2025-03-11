function Figure04_preferenceMapsAcrossAllDatasets(data, fPlots, sets, ...
    retinotopyRF)

limits = [-132 -84 -16 38];
for s = 1:2
    if retinotopyRF(s)
        str = 'retinoptopic';
    else
        str = 'measured';
    end
    % direction
    colors = colmaps.colorcet('C7'); % has 256 entries/colours
    % set directions to range 1-360
    angles = round(mod(data(s).dirPref,360));
    angles(angles == 0) = 360;
    figure
    ind = ~any(isnan(data(s).rfPos),2) & ~isnan(data(s).dirPref);
    scatter(data(s).rfPos(ind,1), data(s).rfPos(ind,2), [], ...
        angles(ind), "filled")
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
    title(sprintf('Preferred directions of %s @ %s RFs (n=%d)', ...
        sets{s}, str, sum(ind)))
    io.saveFigure(gcf, fPlots, sprintf('globalScatter_%s_direction_%sRFs', ...
        sets{s}, str))

    % orientation
    colors = colmaps.colorcet('C1');
    % set orientations to range 1-180
    angles = round(mod(data(s).oriPref,180));
    angles(angles == 0) = 180;
    figure
    ind = ~any(isnan(data(s).rfPos),2) & ~isnan(data(s).oriPref);
    scatter(data(s).rfPos(ind,1), data(s).rfPos(ind,2), [], ...
        angles(ind), "filled")
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
    title(sprintf('Preferred orientations of %s @ %s RFs (n=%d)', ...
        sets{s}, str, sum(ind)))
    io.saveFigure(gcf, fPlots, sprintf('globalScatter_%s_orientation_%sRFs', ...
        sets{s}, str))
end