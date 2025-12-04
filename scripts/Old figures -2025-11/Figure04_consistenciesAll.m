function Figure04_consistenciesAll(maps, fPlots, sets, retinotopyRF, ...
    measures)

%% Parameters
% create surrogate global maps of direction/orientation tuning
numPerm = 1000;
% plotting histograms
binEdges = 0:0.05:1;
bins = binEdges(2:end) - 0.025;
binsSmooth = linspace(bins(1), bins(end), 100);
yLimH = 0.061;

%% Plot histograms + scatters: consistencies compared to null distribution
for s = 1:2
    fprintf('Consistency of %s:\n', sets{s})
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
        xlabel(sprintf('%s consistency', measures{m}))
        ylabel('Probability')
        title(sprintf('%s (%s RFs, n = %d)', sets{s}, str, sum(n)))
        io.saveFigure(gcf, fPlots, sprintf('consistency_%s_%s_histogram_%sRFs', ...
            sets{s}, measures{m}, str))

        mn = mean(maps(s).(measures{m}).consistencies(:), "omitnan");
        sem = std(maps(s).(measures{m}).consistencies(:), "omitnan") / ...
            sum(~isnan(maps(s).(measures{m}).consistencies(:)));
        mnNull = squeeze(mean(maps(s).(measures{m}).nullCons, [1 2], "omitnan"));
        fprintf('  %s:\n', measures{m})
        fprintf('    Original (mean +- SEM): %.4f +- %.4f\n', mn, sem)
        fprintf('    Null (median, 2.5-97.5 interval): %.4f (%.4f-%.4f)\n', ...
            prctile(mnNull, [50 2.5 97.5]))
    end
end