function Figure07_polarHistograms(glob, fPlots, data)

warning('off', 'MATLAB:print:ContentTypeImageSuggested')

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
gridX = -90 : gridSpace : 0;
gridY = 40 : -gridSpace : -10;
polarEdges = deg2rad(-7.5:15:352.5);
cols = {'k', 'b', 'r', 'm'};
% matching vectors
minCount = 4;
% prediction errors
edges_err = 0:5:60;
yLimits = [30 45 75];
heights = [28 43 70];

% Tests
nPerm = 1000;

% SC layers
SCLabel = {'sSC', 'dSC', 'SC'};

%% Optimize pitch angle for each session
% minimize distance between preferred direction/orientation and nearest
% predicted vector
bestPitches = NaN(length(data), 1);
animals = unique({data.animal});
options = optimset('Display', 'off'); % optimset('Display', 'iter');
for subj = 1:length(animals)
    sessions = strcmp(animals{subj}, {data.animal});
    d = cat(1, data(sessions).dirPreferences);
    d(~cat(1, data(sessions).dirTuned)) = NaN;
    o = cat(1, data(sessions).oriPreferences);
    o(~cat(1, data(sessions).oriTuned)) = NaN;
    fun = @(x)algebra.determinePitchError(x, d, o, ...
        cat(1, data(sessions).rfPos));
    pitch = fminbnd(fun, pitchRange(1), pitchRange(2), options);
    bestPitches(sessions) = pitch;
end

%% Pool data across sessions
rfPositions = cat(1, data.rfPos);
directionPreferences = cat(1, data.dirPreferences);
directionPreferences(~cat(1, data.dirTuned)) = NaN;
orientationPreferences = cat(1, data.oriPreferences);
orientationPreferences(~cat(1, data.oriTuned)) = NaN;
ID = [];
layer = [];
modelDirs = []; % in degree
modelOris_long = []; % in degree
modelOris_lat = []; % in degree
for session = 1:length(data)
    n = length(data(session).depth);

    ID = [ID; ones(n,1) .* session];
    l = ones(n, 1);
    l(data(session).depth > data(session).SO_depth) = 2;
    layer = [layer; l];

    tv = NaN(n, 4);
    lov = NaN(n, 2);
    lav = NaN(n, 2);
    [ds_trans, os_long, os_lat] = algebra.getDsOsAxes(bestPitches(session));
    for unit = 1:n
        if any(isnan(data(session).rfPos(unit,:)))
            continue
        end
        if data(session).dirTuned(unit)
            for v = 1:4
                [~, ~, tv(unit,v)] = algebra.getTranslationDir( ...
                    ds_trans(v,:), data(session).rfPos(unit,:));
            end
        end
        if data(session).oriTuned(unit)
            for v = 1:2
                [~, ~, lov(unit,v)] = algebra.getTranslationDir( ...
                    os_long(v,:), data(session).rfPos(unit,:));
                [~, ~, lav(unit,v)] = algebra.getLatitudeOrientation( ...
                    os_lat(v,:), data(session).rfPos(unit,:));
            end
        end
    end
    modelDirs = [modelDirs; tv];
    modelOris_long = [modelOris_long; mod(lov, 180)];
    modelOris_lat = [modelOris_lat; lav];
end

%% Polar histograms
pitch = mean(bestPitches);
[ds_trans, os_long, os_lat] = algebra.getDsOsAxes(pitch);
for l = 3 % layer in SC: (1) sSC, (2) dSC, (3) pooled
    if l < 3
        valid = layer == l;
    else
        valid = true(size(layer));
    end

    % direction
    figure('Position', [1 49 1920 955])
    tiledlayout(length(gridY), length(gridX), "TileSpacing", "tight")
    for x = 1:length(gridX)
        ind0 = rfPositions(:,1) > gridX(x)-0.5*gridSpace & ...
            rfPositions(:,1) <= gridX(x)+0.5*gridSpace;
        for y = 1:length(gridY)
            ind1 = ind0 & rfPositions(:,2) <= gridY(y)+0.5*gridSpace & ...
                rfPositions(:,2) > gridY(y)-0.5*gridSpace;
            dp = deg2rad(directionPreferences(ind1 & valid));
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
    sgtitle(sprintf(['Direction tuning with pitch %.2f deg\n\n' ...
        '[%d to %d azim, %d to %d elev]'], ...
        pitch, gridX(1), gridX(end), gridY(end), gridY(1)))
    io.saveFigure(gcf, fPlots, ...
        sprintf('polarHists_direction_%s', SCLabel{l}))

    % orientation
    figure('Position', [1 49 1920 955])
    tiledlayout(length(gridY), length(gridX), "TileSpacing", "tight")
    for x = 1:length(gridX)
        ind0 = rfPositions(:,1) > gridX(x)-0.5*gridSpace & ...
            rfPositions(:,1) <= gridX(x)+0.5*gridSpace;
        for y = 1:length(gridY)
            ind1 = ind0 & rfPositions(:,2) <= gridY(y)+0.5*gridSpace & ...
                rfPositions(:,2) > gridY(y)-0.5*gridSpace;
            op = deg2rad(orientationPreferences(ind1 & valid));
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
    sgtitle(sprintf(['Orientation tuning with pitch %.2fdeg\n' ...
        'V_{lat} axis: [%d %d]\n[%d to %d azim, %d to %d elev]'], ...
        pitch, round(os_lat(2,:)), gridX([1 end]), gridY([end 1])))
    io.saveFigure(gcf, fPlots, ...
        sprintf('polarHists_orientation_%s', SCLabel{l}))
end

%% Prediction error, significance tests and classification histograms
% prediction error for original data
errDir_original = cell(3,1);
errOri_original = cell(3,1);
for l = 1:2
    valid = ~any(isnan(rfPositions), 2);
    if l < 3
        valid = valid & layer == l;
    end

    edges = .5:4.5;
    bins = 1:4;

    % 1. direction
    ind = valid & ~isnan(directionPreferences);
    dp_sets = ID(ind);
    err = abs(directionPreferences(ind) - modelDirs(ind,:));
    ind = err > 180;
    err(ind) = 360 - err(ind);
    [err, optTransVects] = min(err, [], 2);
    errDir_original{l} = err;

    % plot histograms of closest vectors
    counts = NaN(max(dp_sets), length(bins));
    for k = unique(dp_sets)'
        ind = dp_sets == k;
        counts(k,:) = histcounts(optTransVects(ind), edges);
    end
    n = sum(counts, 2);
    counts(n < minCount,:) = [];
    n(n < minCount) = [];
    figure('Position', glob.figPositionDefault)
    swarmchart(repmat(bins, length(n), 1), counts ./ n, 2 * n, ...
        'filled', 'XJitter', 'density', 'XJitterWidth', 0.5);
    hold on
    plot([-0.4; 0.4] + (1:4), ...
        repmat(median(counts ./ n, 1, "omitnan"), 2, 1), 'k', ...
        'LineWidth', 4)
    ylim([0 1])
    set(gca, 'XTick', 1:4, 'XTickLabel', ...
        {'advance','retreat','rise','fall'}, 'Box', 'off')
    title(sprintf('Direction vectors - %s (n = %d)', SCLabel{l}, sum(~isnan(n))))
    io.saveFigure(gcf, fPlots, sprintf('optimalVectors_direction_%s', SCLabel{l}))
    
    % 2. orientation
    ind = valid & ~isnan(orientationPreferences);
    op_sets = ID(ind);
    err = abs(orientationPreferences(ind) - ...
        [modelOris_long(ind,:), modelOris_lat(ind,:)]);
    ind = err > 90;
    err(ind) = 180 - err(ind);
    [err, optOriVects] = min(err, [], 2);
    errOri_original{l} = err;

    % plot histograms of closest vectors
    counts = NaN(max(op_sets), length(bins));
    for k = unique(op_sets)'
        ind = op_sets == k;
        counts(k,:) = histcounts(optOriVects(ind), edges);
    end
    n = sum(counts, 2);
    counts(n < minCount,:) = [];
    n(n < minCount) = [];
    figure('Position', glob.figPositionDefault)
    swarmchart(repmat(bins, length(n), 1), counts ./ n, 2 * n, ...
        'filled', 'XJitter', 'density', 'XJitterWidth', 0.5);
    hold on
    plot([-0.4; 0.4] + (1:4), ...
        repmat(median(counts ./ n, 1, "omitnan"), 2, 1), 'k', ...
        'LineWidth', 4)
    set(gca, 'XTick', 1:4, 'XTickLabel', ...
        {'H_{long}','V_{long}','H_{lat}','V_{lat}'}, 'Box', 'off')
    ylim([0 1])
    title(sprintf('Orientation vectors - %s (n = %d)', SCLabel{l}, sum(~isnan(n))))
    io.saveFigure(gcf, fPlots, sprintf('optimalVectors_orientation_%s', SCLabel{l}))
end

% prediciton error for permuted data
errDir_perm = NaN(3, nPerm);
errOri_perm = NaN(3, nPerm);
rng('default');
for l = 1:2
    valid = ~any(isnan(rfPositions), 2);
    if l < 3
        valid = valid & layer == l;
    end
    for p = 1:nPerm
        % 1. direction
        ind = find(valid & ~isnan(directionPreferences));
        order = randperm(length(ind));
        dp_perm = directionPreferences(ind(order));
        err_tmp = abs(dp_perm - modelDirs(ind,:));
        ind = err_tmp > 180;
        err_tmp(ind) = 360 - err_tmp(ind);
        err_tmp = min(err_tmp, [], 2);
        errDir_perm(l,p) = median(err_tmp);
        % 2. orientation
        ind = find(valid & ~isnan(orientationPreferences));
        order = randperm(length(ind));
        op_perm = orientationPreferences(ind(order));
        err_tmp = abs(op_perm - [modelOris_long(ind,:), modelOris_lat(ind,:)]);
        ind = err_tmp > 90;
        err_tmp(ind) = 180 - err_tmp(ind);
        err_tmp = min(err_tmp, [], 2);
        errOri_perm(l,p) = median(err_tmp);
    end
end
errDir_perm = prctile(errDir_perm, [2.5 50 97.5], 2);
errOri_perm = prctile(errOri_perm, [2.5 50 97.5], 2);

% prediction error for uniform data
errDir_uni = NaN(3, 3);
errOri_uni = NaN(3, 3);
rng('default');
% 1. direction
valid = find(~isnan(directionPreferences) & ~any(isnan(rfPositions), 2));
dp_uni = rand(length(valid), 1, nPerm) .* 360;
err_tmp = abs(dp_uni - modelDirs(valid,:));
ind = err_tmp > 180;
err_tmp(ind) = 360 - err_tmp(ind);
err_tmp = squeeze(min(err_tmp, [], 2));
errDir_uni(1,:) = prctile( median( err_tmp(layer(valid)==1,:), 1 ), ...
    [2.5 50 97.5]); % sSC
errDir_uni(2,:) = prctile( median( err_tmp(layer(valid)==2,:), 1 ), ...
    [2.5 50 97.5]); % dSC
errDir_uni(3,:) = prctile( median( err_tmp, 1 ), ...
    [2.5 50 97.5]); % SC
% 2. orientation
valid = find(~isnan(orientationPreferences) & ~any(isnan(rfPositions), 2));
op_uni = rand(length(valid), 1, nPerm) .* 180;
err_tmp = abs(op_uni - [modelOris_long(valid,:), modelOris_lat(valid,:)]);
ind = err_tmp > 90;
err_tmp(ind) = 180 - err_tmp(ind);
err_tmp = min(err_tmp, [], 2);
errOri_uni(1,:) = prctile( median( err_tmp(layer(valid)==1,:), 1 ), ...
    [2.5 50 97.5]); % sSC
errOri_uni(2,:) = prctile( median( err_tmp(layer(valid)==2,:), 1 ), ...
    [2.5 50 97.5]); % dSC
errOri_uni(3,:) = prctile( median( err_tmp, 1 ), ...
    [2.5 50 97.5]); % SC

% plot prediction errors of original data and confidence intervals for
% permuted and uniform data
for l = 1:2
    % 1. direction
    h = [0 0];
    figure('Position', glob.figPositionDefault)
    histogram(errDir_original{l}, edges_err, "FaceColor", "k", "FaceAlpha", 1)
    hold on
    plot(median(errDir_original{l}), heights(l), 'v', ...
        "MarkerFaceColor", 'k', "MarkerEdgeColor", "none")
    h(1) = plot(errDir_perm(l, [1 3]), [1 1] .* heights(l), ...
        "Color", [1 1 1] .* 0.4, "LineWidth", 2);
    plot(errDir_perm(l, 2), heights(l), '.', ...
        "Color", [1 1 1] .* 0.4, "MarkerSize", 15)
    h(2) = plot(errDir_uni(l, [1 3]), [1 1] .* heights(l), ...
        "Color", [1 1 1] .* 0.8, "LineWidth", 2);
    plot(errDir_uni(l, 2), heights(l), '.', ...
        "Color", [1 1 1] .* 0.8, "MarkerSize", 15)
    ylim([0 yLimits(l)])
    legend(h, 'permuted', 'uniform')
    set(gca, "Box", "off")
    xlabel('Prediction error (deg)')
    ylabel('#Units')
    title('Preferred directions')
    io.saveFigure(gcf, fPlots, sprintf('predictionErrors_direction_%s', ...
        SCLabel{l}))
    % 2. orientation
    h = [0 0];
    figure('Position', glob.figPositionDefault)
    histogram(errOri_original{l}, edges_err, "FaceColor", "k", "FaceAlpha", 1)
    hold on
    plot(median(errOri_original{l}), heights(l), 'v', ...
        "MarkerFaceColor", 'k', "MarkerEdgeColor", "none")
    h(1) = plot(errOri_perm(l, [1 3]), [1 1] .* heights(l), ...
        "Color", [1 1 1] .* 0.4, "LineWidth", 2);
    plot(errOri_perm(l, 2), heights(l), '.', ...
        "Color", [1 1 1] .* 0.4, "MarkerSize", 15)
    h(2) = plot(errOri_uni(l, [1 3]), [1 1] .* heights(l), ...
        "Color", [1 1 1] .* 0.8, "LineWidth", 2);
    plot(errOri_uni(l, 2), heights(l), '.', ...
        "Color", [1 1 1] .* 0.8, "MarkerSize", 15)
    ylim([0 yLimits(l)])
    legend(h, 'permuted', 'uniform')
    set(gca, "Box", "off")
    xlabel('Prediction error (deg)')
    ylabel('#Units')
    title('Preferred orientations')
    io.saveFigure(gcf, fPlots, sprintf('predictionErrors_orientation_%s', ...
        SCLabel{l}))
end

fprintf('Prediction errors [95%% conf int: (1) permuted data, (2) uniform distribution]:\n')
fprintf('  Direction\n')
for l = 1:2
    fprintf('    %s\n', SCLabel{l})
    fprintf('      Median: %.2f [%.2f-%.2f, %.2f-%.2f]\n', ...
        median(errDir_original{l}), ...
        errDir_perm(l,[1 3]), ...
        errDir_uni(l,[1 3]))
end
fprintf('  Orientation\n')
for l = 1:2
    fprintf('    %s\n', SCLabel{l})
    fprintf('      Median: %.2f [%.2f-%.2f, %.2f-%.2f]\n', ...
        median(errOri_original{l}), ...
        errOri_perm(l,[1 3]), ...
        errOri_uni(l,[1 3]))
end
