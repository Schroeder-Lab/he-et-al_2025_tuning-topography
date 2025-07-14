function Figure05_polarHists(data, fPlots, sets, selectivityThresholds, suffix)

if nargin < 4 
    selectivityThresholds = [0 1; 0 1];
end
if nargin  < 5
    suffix = '';
end

%% Parameters
[ds_trans, os_long, os_lat] = algebra.getDsOsAxes();

% Plotting
gridSpace = 10;
gridX = -130 : gridSpace : -90;
gridY = 40 : -gridSpace : -10;

polarEdges = deg2rad(-5:10:355);
cols = {'k', 'b', 'r', 'm'};

% Tests
nPerm = 1000;

%% Polar histograms
for s = 1:length(data)
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
    sgtitle(sprintf('%s: direction tuning [%d to %d azim, %d to %d elev]', sets{s}, ...
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
    sgtitle(sprintf('%s: orientation tuning [%d to %d azim, %d to %d elev]', sets{s}, ...
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
    transVectors = NaN(length(valid), 4); % 4 direction vectors of DSGCs at RF position
    for k = 1:length(valid)
        for v = 1:4
            [~,~, predDir] = algebra.getTranslationDir(ds_trans(v,:), ...
                data(s).rfPos(valid(k),:));
            transVectors(k,v) = predDir;
        end
    end
    % orientation
    % only use units with minimum OSI and maximum DSI
    valid = find(data(s).OSI >= selectivityThresholds(2,1) & ...
        (data(s).DSI <= selectivityThresholds(2,2) | isnan(data(s).DSI)) & ...
        ~any(isnan(data(s).rfPos), 2));
    op = data(s).oriPref(valid); % preferred orientations
    % tranform angles from direction of motion (axial motion) to
    % orientation of grating
    op = mod(op + 90, 180);
    longVectors = NaN(length(valid), 2); % 2 longitudinal vectors of OSGCs at RF position
    latVectors = NaN(length(valid), 2); % 2 latitudinal vectors of OSGCs at RF position
    for k = 1:length(valid)
        for v = 1:2
            [~,~, predOri] = algebra.getTranslationDir(os_long(v,:), ...
                data(s).rfPos(valid(k),:));
            longVectors(k,v) = mod(predOri, 180);
            [~,~, predOri] = algebra.getLatitudeOrientation(os_lat(v,:), ...
                data(s).rfPos(valid(k),:));
            latVectors(k,v) = predOri;
        end
    end

    % prediction error for original data
    % 1. direction
    err = abs(dp - transVectors);
    ind = err > 180;
    err(ind) = 360 - err(ind);
    [errDir_original{s}, optTransVects] = min(err, [], 2);
    figure
    histogram(optTransVects, .5:4.5)
    set(gca, 'XTick', 1:4, 'XTickLabel', ...
        {'advance','retreat','rise','fall'}, 'Box', 'off')
    title(sprintf('%s: Direction vectors (n = %d)', sets{s}, length(optTransVects)))
    io.saveFigure(gcf, fPlots, sprintf('optimalVectors_%s_direction%s', ...
        sets{s}, suffix))
    % 2. orientation
    err = abs(op - [longVectors, latVectors]);
    ind = err > 90;
    err(ind) = 180 - err(ind);
    [errOri_original{s}, optOriVects] = min(err, [], 2);
    figure
    histogram(optOriVects, .5:4.5)
    set(gca, 'XTick', 1:4, 'XTickLabel', ...
        {'H_{long}','V_{long}','H_{lat}','V_{lat}'}, 'Box', 'off')
    title(sprintf('%s: Orientation vectors (n = %d)', sets{s}, length(optOriVects)))
    io.saveFigure(gcf, fPlots, sprintf('optimalVectors_%s_orientation%s', ...
        sets{s}, suffix))

    % prediciton error for permuted data
    errDir_perm = NaN(length(dp), nPerm);
    errOri_perm = NaN(length(op), nPerm);
    rng('default');
    for p = 1:nPerm
        % 1. direction
        order = randperm(length(dp));
        dp_perm = dp(order);
        err_tmp = abs(dp_perm - transVectors);
        ind = err_tmp > 180;
        err_tmp(ind) = 360 - err_tmp(ind);
        errDir_perm(:,p) = min(err_tmp, [], 2);
        % 2. orientation
        order = randperm(length(op));
        op_perm = op(order);
        err_tmp = abs(op_perm - [longVectors, latVectors]);
        ind = err_tmp > 90;
        err_tmp(ind) = 180 - err_tmp(ind);
        errOri_perm(:,p) = min(err_tmp, [], 2);
    end
    
    % prediciton error for uniform data
    rng('default');
    % 1. direction
    dp_uni = rand(length(dp), 1, nPerm) .* 360;
    err_tmp = abs(dp_uni - transVectors);
    ind = err_tmp > 180;
    err_tmp(ind) = 360 - err_tmp(ind);
    errDir_uni = squeeze(min(err_tmp, [], 2));
    % 2. orientation
    op_uni = rand(length(op), 1, nPerm) .* 180;
    err_tmp = abs(op_uni - [longVectors, latVectors]);
    ind = err_tmp > 90;
    err_tmp(ind) = 180 - err_tmp(ind);
    errOri_uni = squeeze(min(err_tmp, [], 2));

    fprintf('Prediction errors [95%% conf int: (1) permuted data, (2) uniform distribution]:\n')
    fprintf('  %s:\n', sets{s})
    fprintf('    Direction\n')
    fprintf('      Median: %.2f [%.2f-%.2f, %.2f-%.2f]\n', ...
        median(errDir_original{s}), ...
        prctile(median(errDir_perm,1),[2.5 97.5]), ...
        prctile(median(errDir_uni,1),[2.5 97.5]))
    fprintf('    Orientation\n')
    fprintf('      Median: %.2f [%.2f-%.2f, %.2f-%.2f]\n', ...
        median(errOri_original{s}), ...
        prctile(median(errOri_perm,1),[2.5 97.5]), ...
        prctile(median(errOri_uni,1),[2.5 97.5]))
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