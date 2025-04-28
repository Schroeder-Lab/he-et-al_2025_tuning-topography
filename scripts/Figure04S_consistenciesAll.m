function Figure04S_consistenciesAll(maps, fPlots, sets, retinotopyRF, ...
    measures)

%% Parameters
% create surrogate global maps of direction/orientation tuning
numPerm = 1000;
% plotting scatterplots
bins2 = 0.01:0.02:1;
yLimS = [200 120];

%% Plot histograms + scatters: consistencies compared to null distribution
for s = 1:2
    if retinotopyRF(s)
        str = 'retinoptopic';
    else
        str = 'measured';
    end
    for m = 1:2
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
        ylim([0 yLimS(s)])
        xlabel('Consistency')
        ylabel('#units per patch')
        title(sprintf('%s consistency (%s RFs) - %s (n = %d)', ...
            measures{m}, str, sets{s}, ...
            sum(~isnan(maps(s).(measures{m}).consistencies), "all")))
        io.saveFigure(gcf, fPlots, sprintf('consistency_%s_%s_scatter-count_%sRFs', ...
            sets{s}, measures{m}, str))
    end
end