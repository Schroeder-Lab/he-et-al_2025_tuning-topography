function Figure05_consistenciesAll(maps, fPlots, sets, retinotopyRF, ...
    measures)

%% Parameters
% create surrogate global maps of direction/orientation tuning
numPerm = 1000;
% plotting histograms
binEdges = 0:0.05:1;
bins = binEdges(2:end) - 0.025;
binsSmooth = linspace(bins(1), bins(end), 100);
yLimH = 0.061;
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
        n = histcounts(maps(s).(measures{m}).consistencies(:), binEdges);
        ns = interp1(bins, n, binsSmooth, "pchip");
        ns = ns ./ sum(ns);
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
        h(2) = plot(binsSmooth, ns, 'r', "LineWidth", 1);
        legend(h, 'null', 'original')
        ylim([0 yLimH])
        xlabel('Consistency')
        ylabel('Probability')
        title(sprintf('%s consistency (%s RFs) - %s (n = %d)', ...
            measures{m}, str, sets{s}, sum(n)))
        io.saveFigure(gcf, fPlots, sprintf('consistency_%s_%s_histogram_%sRFs', ...
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