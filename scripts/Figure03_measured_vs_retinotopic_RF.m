function Figure03_measured_vs_retinotopic_RF(folders, sets, fPlots)

%% Parameters
% for evaluation of receptive fields (significance/goodness)
minEV = 0.01;
minPeak = 5;
% plotting
minUnits = 10;
maxDist = 20;
maxDiam = 30;

%% Plot histogram of distances
for s = 1:2 % boutons and neurons
    dist_rf_pred = [];
    diam_rf = []; % = 2 STDs (STD is mean across x and y)

    subjDirs = dir(fullfile(folders.data, sets{s}, 'SS*'));
    for subj = 1:length(subjDirs) % animals
        name = subjDirs(subj).name;
        dateDirs = dir(fullfile(folders.data, sets{s}, name, '2*'));
        for dt = 1:length(dateDirs) %dates
            date = dateDirs(dt).name;
            f = fullfile(folders.data, sets{s}, name, date);
            % ignore session if visual noise stimuolus was not present
            if ~isfile(fullfile(f, '_ss_rf.posRetinotopy.npy'))
                continue
            end

            % load data
            rfData = io.getRFFits(f);
            fitPars = rfData.fitParameters;
            rfPos = fitPars(:, [2 4]); % (azimuth, elevation) in visual degrees
            ev_rf = rfData.EV;
            rf_peaks = rfData.peaks;
            outliers = readNPY(fullfile(f, '_ss_rf.outliers.npy'));
            pred = readNPY(fullfile(f, '_ss_rf.posRetinotopy.npy'));

            % only consider significant RFs
            valid = ev_rf >= minEV & rf_peaks >= minPeak;
            valid = valid & ~outliers;
            if sum(valid) < minUnits
                continue
            end

            % Eucleadian distances
            d = sqrt( ...
                sum((rfPos(valid,:) - pred(valid,:)) .^ 2, 2) );
            dist_rf_pred = [dist_rf_pred; d];

            % RF diameters (2 mean STDs)
            sz = 2 * mean(fitPars(valid,[3 5]), 2);
            diam_rf = [diam_rf; sz];
        end
    end

    m = mean(dist_rf_pred);
    dist_rf_pred(dist_rf_pred > maxDist) = maxDist;
    figure
    h = histogram(dist_rf_pred, 0:0.5:maxDist);
    h.FaceColor = [1 1 1] .* 0.5;
    hold on
    plot(m, 5, 'v')
    xlabel('Distance: mapped vs. retinotopic (visual degrees)')
    ylabel(sprintf('%s', sets{s}))
    xlim([-.2 maxDist+.2])
    set(gca, "Box", "off")
    title(sprintf('mean: %.2f, n = %d', m, length(dist_rf_pred)))
    io.saveFigure(gcf, fPlots, ...
        sprintf('rf_%s_distance-mapped_retinotopic', sets{s}))

    m = mean(diam_rf);
    diam_rf(diam_rf > maxDiam) = maxDiam;
    figure
    h = histogram(diam_rf, 0:0.5:maxDiam);
    h.FaceColor = [1 1 1] .* 0.5;
    hold on
    plot(m, 5, 'v')
    xlabel('RF diameter (visual degrees)')
    ylabel(sprintf('%s', sets{s}))
    xlim([-.2 maxDiam+.2])
    set(gca, "Box", "off")
    title(sprintf('mean: %.2f, n = %d', m, length(dist_rf_pred)))
    io.saveFigure(gcf, fPlots, ...
        sprintf('rf_%s_diameter', sets{s}))
end