function Figure03_polarHists(glob, fPlots, data, sets, ...
    selectivityThresholds, suffix)

warning('off', 'MATLAB:print:ContentTypeImageSuggested')

if nargin < 5 
    selectivityThresholds = [0 1; 0 1];
end
if nargin  < 6
    suffix = '';
end

%% Parameters
% to optimize pitch angle;
% translational and latitude axes determined by Sabbah et al assume that
% the lambda-bregma axis of the skull is at an angle of 29 deg (bregma
% lower); the angle of our mice is typically somewhat above 0; to match the
% two angles, we will need to tilt the translational and latitude axes
% upwards by 0-29 deg
pitchRange = [-29 0];

% Plotting
gridSpace = 10;
gridX = -130 : gridSpace : -90;
gridY = 40 : -gridSpace : -10;
polarEdges = deg2rad(-5:10:355);
cols = {'k', 'b', 'r', 'm'};
% matching vectors
minCount = 20;
% prediction errors
edges_err = 0:5:60;
yLimit = 350;
height = 250;

% Tests
nPerm = 1000;

%% Optimize pitch angle for each session
% minimize distance between preferred direction/orientation and nearest
% predicted vector
bestPitches = cell(1,2);
modelDirs = cell(1,2);
modelOris = cell(1,2);
for s = 1:length(data)
    animals = unique(data(s).animal);
    bestPitches{s} = NaN(length(animals), 1);
    modelDirs{s} = NaN(length(data(s).dirPref), 4);
    modelOris{s} = NaN(length(data(s).dirPref), 4);
    options = optimset('Display', 'off'); % optimset('Display', 'iter');
    for subj = 1:length(animals)
        units = find(strcmp(animals{subj}, data(s).animal) & ...
            ~any(isnan(data(s).rfPos), 2));
        d = data(s).dirPref(units);
        o = data(s).oriPref(units);
        fun = @(x)algebra.determinePitchError(x, d, o, ...
            data(s).rfPos(units,:));
        pitch = fminbnd(fun, pitchRange(1), pitchRange(2), options);
        bestPitches{s}(subj) = pitch;

        [ds_trans, os_long, os_lat] = algebra.getDsOsAxes(pitch);
        for k = 1:length(units)
            unit = units(k);
            if ~isnan(data(s).dirPref(unit))
                d = NaN(1,4);
                for v = 1:4
                    [~, ~, d(v)] = algebra.getTranslationDir( ...
                        ds_trans(v,:), data(s).rfPos(unit,:));
                end
                modelDirs{s}(unit,:) = d;
            end
            if ~isnan(data(s).oriPref(unit))
                o1 = NaN(1,2);
                o2 = NaN(1,2);
                for v = 1:2
                    [~, ~, o1(v)] = algebra.getTranslationDir( ...
                        os_long(v,:), data(s).rfPos(unit,:));
                    [~, ~, o2(v)] = algebra.getLatitudeOrientation( ...
                        os_lat(v,:), data(s).rfPos(unit,:));
                end
                modelOris{s}(unit,:) = [mod(o1, 180) o2];
            end
        end
    end
end

%% Polar histograms
for s = 1:length(data)
    % get average best pitch
    avgPitch = mean(bestPitches{s});
    [ds_trans, os_long, os_lat] = algebra.getDsOsAxes(avgPitch);

    % direction
    % only use units with minimum DSI and maximum OSI
    valid = data(s).DSI >= selectivityThresholds(1,1) & ...
        (data(s).OSI <= selectivityThresholds(1,2) | isnan(data(s).OSI));
    figure('Position', [1 49 1920 955])
    tiledlayout(length(gridY), length(gridX), "TileSpacing", "tight")
    for x = 1:length(gridX)
        ind0 = data(s).rfPos(:,1) > gridX(x)-0.5*gridSpace & ...
            data(s).rfPos(:,1) <= gridX(x)+0.5*gridSpace;
        for y = 1:length(gridY)
            ind1 = ind0 & data(s).rfPos(:,2) <= gridY(y)+0.5*gridSpace & ...
                data(s).rfPos(:,2) > gridY(y)-0.5*gridSpace;
            dp = deg2rad(data(s).dirPref(ind1 & valid));
            nexttile((y-1)*length(gridX) + x)
            polarhistogram(dp, polarEdges, 'Normalization', ...
                'count', 'FaceColor', 'k', 'EdgeColor', 'none');
            ln = get(gca, "RLim");
            ln = ln(2);
            hold on
            for d = 1:4
                [~,theta] = algebra.getTranslationDir(ds_trans(d,:), ...
                    [gridX(x) gridY(y)]);
                polarplot([0 theta], [0 ln], cols{d}, 'LineWidth', 1)
            end
            title(sprintf('%d', sum(~isnan(dp))))
            set(gca,'ThetaDir','clockwise', 'RTick', [], 'ThetaTickLabel', {})
        end
    end
    sgtitle(sprintf(['%s: direction tuning with pitch %.2f deg\n\n' ...
        '[%d to %d azim, %d to %d elev]'], ...
        sets{s}, avgPitch, ...
        gridX(1), gridX(end), gridY(end), gridY(1)))
    io.saveFigure(gcf, fPlots, sprintf('polarHists_%s_direction%s', sets{s}, suffix))

    % orientation
    % only use units with minimum OSI and maximum DSI
    valid = data(s).OSI >= selectivityThresholds(2,1) & ...
        (data(s).DSI <= selectivityThresholds(2,2) | isnan(data(s).DSI));
    figure('Position', [1 49 1920 955])
    tiledlayout(length(gridY), length(gridX), "TileSpacing", "tight")
    for x = 1:length(gridX)
        ind0 = data(s).rfPos(:,1) > gridX(x)-0.5*gridSpace & ...
            data(s).rfPos(:,1) <= gridX(x)+0.5*gridSpace;
        for y = 1:length(gridY)
            ind1 = ind0 & data(s).rfPos(:,2) <= gridY(y)+0.5*gridSpace & ...
                data(s).rfPos(:,2) > gridY(y)-0.5*gridSpace;
            op = deg2rad(data(s).oriPref(ind1 & valid));
            % duplicate the histogram for the opposite hemicircle
            op = [op; op + pi];
            % tranform angles from direction of motion (axial motion) to
            % orientation of grating so that bars in histograms align with
            % orientation of grating
            op = op + pi/2;
            nexttile((y-1)*length(gridX) + x)
            polarhistogram(op, polarEdges, 'Normalization', ...
                'count', 'FaceColor', 'k', 'EdgeColor', 'none');
            ln = get(gca, "RLim");
            ln = ln(2);
            hold on
            % longitudes
            for d = 1:2
                [~,theta] = algebra.getTranslationDir(os_long(d,:), ...
                    [gridX(x) gridY(y)]);
                polarplot([theta theta + pi], [ln ln], cols{d}, 'LineWidth', 1)
            end
            % latitudes
            for d = 1:2
                [~,theta] = algebra.getLatitudeOrientation(os_lat(d,:), ...
                    [gridX(x) gridY(y)]);
                polarplot([theta theta + pi], [ln ln], cols{d+2}, 'LineWidth', 1)
            end
            title(sprintf('%d', sum(~isnan(op))/2))
            set(gca,'ThetaDir','clockwise', 'RTick', [], 'ThetaTickLabel', {})
        end
    end
    sgtitle(sprintf(['%s: orientation tuning with pitch %.2fdeg\n' ...
        'V_{lat} axis: [%d %d]\n[%d to %d azim, %d to %d elev]'], ...
        sets{s}, avgPitch, round(os_lat(2,:)), ...
        gridX(1), gridX(end), gridY(end), gridY(1)))
    io.saveFigure(gcf, fPlots, sprintf('polarHists_%s_orientation%s', sets{s}, suffix))
end

%% Prediction error, significance tests and classification histograms
errDir_original = cell(1,2);
errOri_original = cell(1,2);
for s = 1:length(data)
    % direction
    % only use units with minimum DSI and maximum OSI
    valid = find(data(s).DSI >= selectivityThresholds(1,1) & ...
        (data(s).OSI <= selectivityThresholds(1,2) | isnan(data(s).OSI)) & ...
        ~any(isnan(data(s).rfPos), 2));
    dp = data(s).dirPref(valid); % preferred directions
    dp_sets = data(s).set(valid); % IDs of datasets
    predictedDirs = modelDirs{s}(valid,:); % 4 direction vectors of DSGCs at RF position
    % orientation
    % only use units with minimum OSI and maximum DSI
    valid = find(data(s).OSI >= selectivityThresholds(2,1) & ...
        (data(s).DSI <= selectivityThresholds(2,2) | isnan(data(s).DSI)) & ...
        ~any(isnan(data(s).rfPos), 2));
    op = data(s).oriPref(valid); % preferred orientations
    op_sets = data(s).set(valid); % IDs of datasets
    % tranform angles from direction of motion (axial motion) to
    % orientation of grating
    op = mod(op + 90, 180);
    predictedOris = modelOris{s}(valid,:);

    % prediction error for original data
    % 1. direction
    err = abs(dp - predictedDirs);
    ind = err > 180;
    err(ind) = 360 - err(ind);
    [errDir_original{s}, optTransVects] = min(err, [], 2);
    % 2. orientation
    err = abs(op - predictedOris);
    ind = err > 90;
    err(ind) = 180 - err(ind);
    [errOri_original{s}, optOriVects] = min(err, [], 2);

    % prediction error for permuted data
    errDir_perm = NaN(1, nPerm);
    errOri_perm = NaN(1, nPerm);
    rng('default');
    for p = 1:nPerm
        % 1. direction
        order = randperm(length(dp));
        dp_perm = dp(order);
        err_tmp = abs(dp_perm - predictedDirs);
        ind = err_tmp > 180;
        err_tmp(ind) = 360 - err_tmp(ind);
        err_tmp = min(err_tmp, [], 2);
        errDir_perm(p) = median(err_tmp);
        % 2. orientation
        order = randperm(length(op));
        op_perm = op(order);
        err_tmp = abs(op_perm - predictedOris);
        ind = err_tmp > 90;
        err_tmp(ind) = 180 - err_tmp(ind);
        err_tmp = min(err_tmp, [], 2);
        errOri_perm(p) = median(err_tmp);
    end
    errDir_perm = prctile(errDir_perm, [2.5 50 97.5]);
    errOri_perm = prctile(errOri_perm, [2.5 50 97.5]);
    
    % prediciton error for uniform data
    rng('default');
    % 1. direction
    dp_uni = rand(length(dp), 1, nPerm) .* 360;
    err_tmp = abs(dp_uni - predictedDirs);
    ind = err_tmp > 180;
    err_tmp(ind) = 360 - err_tmp(ind);
    err_tmp = squeeze(min(err_tmp, [], 2));
    errDir_uni = prctile(median(err_tmp, 1), [2.5 50 97.5]);
    % 2. orientation
    op_uni = rand(length(op), 1, nPerm) .* 180;
    err_tmp = abs(op_uni - predictedOris);
    ind = err_tmp > 90;
    err_tmp(ind) = 180 - err_tmp(ind);
    err_tmp = squeeze(min(err_tmp, [], 2));
    errOri_uni = prctile(median(err_tmp, 1), [2.5 50 97.5]);

    % Plot histograms of closest vectors
    edges = .5:4.5;
    bins = 1:4;
    % 1. direction
    counts = NaN(max(dp_sets), length(bins));
    for k = unique(dp_sets)'
        ind = dp_sets == k;
        counts(k,:) = histcounts(optTransVects(ind), edges);
    end
    n = sum(counts, 2);
    counts(n < minCount,:) = [];
    n(n < minCount) = [];
    figure('Position', glob.figPositionDefault)
    swarmchart(repmat(bins, length(n), 1), counts ./ n, n, ...
        'filled', 'XJitter', 'density', 'XJitterWidth', 0.5);
    hold on
    plot([-0.4; 0.4] + (1:4), ...
        repmat(median(counts ./ n, 1, "omitnan"), 2, 1), 'k', ...
        'LineWidth', 4)
    ylim([0 1])
    set(gca, 'XTick', bins, 'XTickLabel', ...
        {'advance','retreat','rise','fall'}, 'Box', 'off')
    title(sprintf('%s: Direction vectors (n = %d)', sets{s}, sum(~isnan(n))))
    io.saveFigure(gcf, fPlots, sprintf('optimalVectors_%s_direction%s', ...
        sets{s}, suffix))
    % 2. orientation
    counts = NaN(max(op_sets), length(bins));
    for k = unique(op_sets)'
        ind = op_sets == k;
        counts(k,:) = histcounts(optOriVects(ind), edges);
    end
    n = sum(counts, 2);
    counts(n < minCount,:) = [];
    n(n < minCount) = [];
    figure('Position', glob.figPositionDefault)
    swarmchart(repmat(bins, length(n), 1), counts ./ n, n, ...
        'filled', 'XJitter', 'density', 'XJitterWidth', 0.5);
    hold on
    plot([-0.4; 0.4] + (1:4), ...
        repmat(median(counts ./ n, 1, "omitnan"), 2, 1), 'k', ...
        'LineWidth', 4)
    ylim([0 1])
    set(gca, 'XTick', bins, 'XTickLabel', ...
        {'H_{long}','V_{long}','H_{lat}','V_{lat}'}, 'Box', 'off')
    title(sprintf('%s: Orientation vectors (n = %d)', sets{s}, sum(~isnan(n))))
    io.saveFigure(gcf, fPlots, sprintf('optimalVectors_%s_orientation%s', ...
        sets{s}, suffix))

    % plot prediction errors of original data and confidence intervals for
    % permuted and uniform data
    % 1. direction
    h = [0 0];
    figure('Position', glob.figPositionDefault)
    histogram(errDir_original{s}, edges_err, "FaceColor", "k", "FaceAlpha", 1)
    hold on
    plot(median(errDir_original{s}), height, 'v', ...
        "MarkerFaceColor", 'k', "MarkerEdgeColor", "none")
    h(1) = plot(errDir_perm([1 3]), [1 1] .* height, ...
        "Color", [1 1 1] .* 0.4, "LineWidth", 2);
    plot(errDir_perm(2), height, '.', ...
        "Color", [1 1 1] .* 0.4, "MarkerSize", 15)
    h(2) = plot(errDir_uni([1 3]), [1 1] .* height, ...
        "Color", [1 1 1] .* 0.8, "LineWidth", 2);
    plot(errDir_uni(2), height, '.', ...
        "Color", [1 1 1] .* 0.8, "MarkerSize", 15)
    ylim([0 yLimit])
    legend(h, 'permuted', 'uniform')
    set(gca, "Box", "off")
    xlabel('Prediction error (deg)')
    ylabel(['#' sets{s}])
    title('Preferred directions')
    io.saveFigure(gcf, fPlots, sprintf('predictionErrors_%s_directions%s', ...
        sets{s}, suffix))
    % 2. orientation
    h = [0 0];
    figure('Position', glob.figPositionDefault)
    histogram(errOri_original{s}, edges_err, "FaceColor", "k", "FaceAlpha", 1)
    hold on
    plot(median(errOri_original{s}), height, 'v', ...
        "MarkerFaceColor", 'k', "MarkerEdgeColor", "none")
    h(1) = plot(errOri_perm([1 3]), [1 1] .* height, ...
        "Color", [1 1 1] .* 0.4, "LineWidth", 2);
    plot(errOri_perm(2), height, '.', ...
        "Color", [1 1 1] .* 0.4, "MarkerSize", 15)
    h(2) = plot(errOri_uni([1 3]), [1 1] .* height, ...
        "Color", [1 1 1] .* 0.8, "LineWidth", 2);
    plot(errOri_uni(2), height, '.', ...
        "Color", [1 1 1] .* 0.8, "MarkerSize", 15)
    ylim([0 yLimit])
    legend(h, 'permuted', 'uniform')
    set(gca, "Box", "off")
    xlabel('Prediction error (deg)')
    ylabel(['#' sets{s}])
    title('Preferred orientations')
    io.saveFigure(gcf, fPlots, sprintf('predictionErrors_%s_orientations%s', ...
        sets{s}, suffix))
    
    fprintf('Prediction errors [95%% conf int: (1) permuted data, (2) uniform distribution]:\n')
    fprintf('  %s:\n', sets{s})
    fprintf('    Direction\n')
    fprintf('      Median: %.2f [%.2f-%.2f, %.2f-%.2f]\n', ...
        median(errDir_original{s}), ...
        errDir_perm([1 3]), ...
        errDir_uni([1 3]))
    fprintf('    Orientation\n')
    fprintf('      Median: %.2f [%.2f-%.2f, %.2f-%.2f]\n', ...
        median(errOri_original{s}), ...
        errOri_perm([1 3]), ...
        errOri_uni([1 3]))
end

fprintf('Prediction error: boutons vs neurons:\n')
fprintf('  Direction:\n')
fprintf('    Boutons - neurons (median): %.2f (p = %.4f, Wilcoxon rank-sum)\n', ...
    median(errDir_original{1})-median(errDir_original{2}), ...
    ranksum(errDir_original{1}, errDir_original{2}))
fprintf('  Orientation:\n')
fprintf('    Boutons - neurons (median): %.2f (p = %.4f, Wilcoxon rank-sum)\n', ...
    median(errOri_original{1})-median(errOri_original{2}), ...
    ranksum(errOri_original{1}, errOri_original{2}))