function Figure02_measured_vs_retinotopic_RF(folders, glob, sets, fPlots)

%% Parameters
% for evaluation of receptive fields (significance/goodness)
minEV = 0.01;
minPeak = 5;
% plotting
minUnits_RF = 10;
maxDist = 20;
maxDiam = 30;
colSets = lines(2);
edgesDist = 0:0.5:maxDist;
edgesSize = 0:0.5:maxDiam;

%% Plot histogram of RF sizes and distances (mapped vs retinotopic)
figDist = figure('Position', glob.figPositionDefault);
hold on
hDist = [0 0];
mDist = [0 0];
numDist = [0 0];
figSize = figure('Position', glob.figPositionDefault);
hold on
hSize = [0 0];
numSize = [0 0];
numTotal = [0 0];
numSessions = [0 0];
numAnimals = [0 0];
dist_rf_pred = cell(2,1);

% from RFs corrected for eye position
mSize_corr = [0 0];
numSize_corr = [0 0];
numTotal_corr = [0 0];
numSessions_corr = [0 0];
diam_rf_corr = cell(2,1); % = 2 STDs (STD is mean across x and y)

for s = 1:2 % boutons and neurons
    subjDirs = dir(fullfile(folders.data, sets{s}, 'SS*'));
    for subj = 1:length(subjDirs) % animals
        name = subjDirs(subj).name;
        dateDirs = dir(fullfile(folders.data, sets{s}, name, '2*'));
        for dt = 1:length(dateDirs) %dates
            date = dateDirs(dt).name;
            f = fullfile(folders.data, sets{s}, name, date);
            % ignore session if visual noise stimulus was not present
            if ~isfile(fullfile(f, '_ss_rf.posRetinotopy.npy'))
                continue
            end

            % load data (RFs not corrected for eye position)
            rfData = io.getNoiseRFFits(f, false);
            fitPars = rfData.gaussPars;
            rfPos = fitPars(:, [2 4]); % (azimuth, elevation) in visual degrees
            ev_rf = rfData.EV;
            rf_peaks = rfData.peaks;
            outliers = readNPY(fullfile(f, '_ss_rf.outliers.npy'));
            pred = readNPY(fullfile(f, '_ss_rf.posRetinotopy.npy'));

            % only consider significant RFs
            valid = ev_rf >= minEV & rf_peaks >= minPeak;
            valid = valid & ~outliers;
            if sum(valid) < minUnits_RF
                continue
            end
            numTotal(s) = numTotal(s) + length(ev_rf);
            numSessions(s) = numSessions(s) + 1;

            % Eucleadian distances
            d = sqrt( ...
                sum((rfPos(valid,:) - pred(valid,:)) .^ 2, 2) );
            dist_rf_pred{s} = [dist_rf_pred{s}; d];


            % load data (RFs corrected for eye position)
            if ~isfile(fullfile(f, '_ss_rf_eye.maps.npy'))
                continue
            end
            rfData = io.getNoiseRFFits(f, true);
            fitPars = rfData.gaussPars;
            ev_rf = rfData.EV;
            rf_peaks = rfData.peaks;

            % only consider significant RFs
            valid = ev_rf >= minEV & rf_peaks >= minPeak;
            valid = valid & ~outliers;
            if sum(valid) < minUnits_RF
                continue
            end
            numTotal_corr(s) = numTotal_corr(s) + length(ev_rf);
            numSessions_corr(s) = numSessions_corr(s) + 1;

            % RF diameters (2 mean STDs)
            sz = 2 * mean(fitPars(valid,[3 5]), 2);
            diam_rf_corr{s} = [diam_rf_corr{s}; sz];
        end
        numAnimals(s) = numAnimals(s) + 1;
    end

    numSize(s) = length(dist_rf_pred{s});
    mDist(s) = median(dist_rf_pred{s});
    numDist(s) = length(dist_rf_pred{s});
    dist_rf_pred{s}(dist_rf_pred{s} > maxDist) = maxDist;
    figure(figDist)
    n = histcounts(dist_rf_pred{s}, edgesDist, "Normalization", "probability");
    hDist(s) = plot(edgesDist(1:end-1)+0.25, n, "Color", colSets(s,:), "LineWidth", 2);
    plot(mDist(s), 0.14, 'v', "MarkerFaceColor", colSets(s,:), ...
        "MarkerEdgeColor", "none")

    numSize_corr(s) = length(diam_rf_corr{s});
    mSize_corr(s) = median(diam_rf_corr{s});
    diam_rf_corr{s}(diam_rf_corr{s} > maxDiam) = maxDiam;
    figure(figSize)
    n = histcounts(diam_rf_corr{s}, edgesSize, "Normalization", "probability");
    hSize(s) = plot(edgesSize(1:end-1)+0.25, n, "Color", colSets(s,:), "LineWidth", 2);
    plot(mSize_corr(s), 0.12, 'v', "MarkerFaceColor", colSets(s,:), ...
        "MarkerEdgeColor", "none")
end

figure(figDist)
xlim([-.2 maxDist+.2])
set(gca, "Box", "off")
l = legend(hDist, ...
    sprintf('%s (median: %.2f, n = %d)', sets{1}, mDist(1), numDist(1)), ...
    sprintf('%s (median: %.2f, n = %d)', sets{2}, mDist(2), numDist(2)));
l.Box = "off";
xlabel('Distance: mapped vs. retinotopic (visual degrees)')
ylabel('Proportion of units')
% io.saveFigure(gcf, fPlots, 'rf_distance-mapped_retinotopic')

figure(figSize)
xlim([-.2 maxDiam+.2])
set(gca, "Box", "off")
l = legend(hSize, ...
    sprintf('%s (median: %.2f, n = %d)', sets{1}, mSize_corr(1), numSize_corr(1)), ...
    sprintf('%s (median: %.2f, n = %d)', sets{2}, mSize_corr(2), numSize_corr(2)));
l.Box = "off";
xlabel('RF diameter (visual degrees)')
ylabel('Proportion of units')
io.saveFigure(gcf, fPlots, 'rf_diameter')

%% Tests
fprintf('Number of mapped RFs (uncorrected) of total:\n')
fprintf('  Boutons: %d of %d (%.1f%%), %d sessions, %d animals\n', ...
    numSize(1), numTotal(1), numSize(1) / numTotal(1) * 100, ...
    numSessions(1), numAnimals(1))
fprintf('  Neurons: %d of %d (%.1f%%), %d sessions, %d animals\n', ...
    numSize(2), numTotal(2), numSize(2) / numTotal(2) * 100, ...
    numSessions(2), numAnimals(2))

p = ranksum(dist_rf_pred{1}, dist_rf_pred{2});
fprintf('Distance: mapped vs retinotopic RFs\n')
fprintf('  Boutons: %.4f median\n', mDist(1))
fprintf('  Neurons: %.4f median\n', mDist(2))
fprintf('  Boutons - neurons (median): %.4f (p = %.6f)\n', ...
    median(dist_rf_pred{1}) - median(dist_rf_pred{2}), p)


fprintf('Number of mapped RFs (corrected for eye position) of total:\n')
fprintf('  Boutons: %d of %d (%.1f%%), %d sessions\n', ...
    numSize_corr(1), numTotal_corr(1), numSize_corr(1) / numTotal_corr(1) * 100, ...
    numSessions_corr(1))
fprintf('  Neurons: %d of %d (%.1f%%), %d sessions\n', ...
    numSize_corr(2), numTotal_corr(2), numSize_corr(2) / numTotal_corr(2) * 100, ...
    numSessions(2))

p = ranksum(diam_rf_corr{1}, diam_rf_corr{2});
fprintf('RF sizes\n')
fprintf('  Boutons: %.4f median\n', mSize_corr(1))
fprintf('  Neurons: %.4f median\n', mSize_corr(2))
fprintf('  Boutons - neurons (median): %.4f (p = %.6f)\n', ...
    median(diam_rf_corr{1}) - median(diam_rf_corr{2}), p)