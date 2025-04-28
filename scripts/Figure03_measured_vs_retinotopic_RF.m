function Figure03_measured_vs_retinotopic_RF(folders, sets, fPlots)

%% Parameters
% for evaluation of receptive fields (significance/goodness)
minEV = 0.01;
minPeak = 5;
% plotting
minUnits = 10;
maxDist = 20;
maxDiam = 30;
colSets = lines(2);
edgesDist = 0:0.5:maxDist;
edgesSize = 0:0.5:maxDiam;

%% Plot histogram of RF sizes and distances (mapped vs retinotopic)
figDist = figure;
hold on
hDist = [0 0];
mDist = [0 0];
numDist = [0 0];
figSize = figure;
hold on
hSize = [0 0];
mSize = [0 0];
numSize = [0 0];
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

    mDist(s) = mean(dist_rf_pred);
    numDist(s) = length(dist_rf_pred);
    dist_rf_pred(dist_rf_pred > maxDist) = maxDist;
    figure(figDist)
    n = histcounts(dist_rf_pred, edgesDist, "Normalization", "probability");
    hDist(s) = plot(edgesDist(1:end-1)+0.25, n, "Color", colSets(s,:), "LineWidth", 2);
    % hSize(s) = histogram(dist_rf_pred, 0:0.5:maxDist, "EdgeColor", colSets(s,:), ...
    %     "Normalization", "probability", "DisplayStyle", "stairs");
    % h.FaceColor = [1 1 1] .* 0.5;
    plot(mDist(s), 0.14, 'v', "MarkerFaceColor", colSets(s,:), ...
        "MarkerEdgeColor", "none")

    mSize(s) = mean(diam_rf);
    numSize(s) = length(diam_rf);
    diam_rf(diam_rf > maxDiam) = maxDiam;
    figure(figSize)
    n = histcounts(diam_rf, edgesSize, "Normalization", "probability");
    hSize(s) = plot(edgesSize(1:end-1)+0.25, n, "Color", colSets(s,:), "LineWidth", 2);
    % hSize(s) = histogram(diam_rf, binsSize, "EdgeColor", colSets(s,:), ...
    %     "Normalization", "probability", "DisplayStyle", "stairs");
    % h.FaceColor = [1 1 1] .* 0.5;
    plot(mSize(s), 0.12, 'v', "MarkerFaceColor", colSets(s,:), ...
        "MarkerEdgeColor", "none")
end

figure(figDist)
xlim([-.2 maxDist+.2])
set(gca, "Box", "off")
l = legend(hDist, ...
    sprintf('%s (mean: %.2f, n = %d)', sets{1}, mDist(1), numDist(1)), ...
    sprintf('%s (mean: %.2f, n = %d)', sets{2}, mDist(2), numDist(2)));
l.Box = "off";
xlabel('Distance: mapped vs. retinotopic (visual degrees)')
ylabel('Proportion of units')
io.saveFigure(gcf, fPlots, 'rf_distance-mapped_retinotopic')

figure(figSize)
xlim([-.2 maxDiam+.2])
set(gca, "Box", "off")
l = legend(hSize, ...
    sprintf('%s (mean: %.2f, n = %d)', sets{1}, mSize(1), numSize(1)), ...
    sprintf('%s (mean: %.2f, n = %d)', sets{2}, mSize(2), numSize(2)));
l.Box = "off";
xlabel('RF diameter (visual degrees)')
ylabel('Proportion of units')
io.saveFigure(gcf, fPlots, 'rf_diameter')