function Figure04_polarHists(data, fPlots, sets, selectivityThresholds, suffix)

if nargin < 4 
    selectivityThresholds = [0 1; 0 1];
end
if nargin  < 5
    suffix = '';
end

%% Parameters
% Email from Shai Sabbah, 03.04.2025, all axes in [azimuth, elevation], 
% all axes provided for right eye!
% For DSGCs:
% N, 0, 15
% T, 180, -11
% D, -113, 83
% V, 74, -88
% For OSGCs:
% Longitudinal system:
% H, -7, 1
% V, 82, -64
% Latitudinal system:
% H, 105, -68
% V, -53, -6 (opposite: [127 6])

% Translation vectors for DSGCs from Sabbah et al (2017)
% [azimuth, elevation]
ds_N = [180 -15]; % advance cells, opposite to [0 15] (black)
ds_T = [0 11]; % retreat cells, opposite to [180 -11] (blue)
ds_D = [-67 -83]; % rise cells, opposite to [113 83] (red)
ds_V = [106 88]; % fall cells, opposite to [-74 -88] (magenta)
ds_trans = [ds_N; ds_T; ds_D; ds_V];

% Axes for OSGCs from Laniado et al (2025)
% longitudinal cells
os_long_H = [7 1]; % black
os_long_V = [-82 -64]; % blue
os_long = [os_long_H; os_long_V];
% latitudinal cells
os_lat_H = [-105 -68]; % red
os_lat_V = [53 -6]; % magenta (opposite: [-127 6])
os_lat = [os_lat_H; os_lat_V];

% Plotting
gridSpace = 10;
gridX = -130 : gridSpace : -90;
gridY = 40 : -gridSpace : -10;

polarEdges = deg2rad(-5:10:355);
cols = {'k', 'b', 'r', 'm'};

% Tests
nPerm = 1000;

%% Make figures
for s = 1:length(data)
    % direction
    % only use units with minimum DSI and maximum OSI
    valid = data(s).DSI >= selectivityThresholds(1,1) & ...
        data(s).OSI <= selectivityThresholds(1,2);
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
        data(s).DSI <= selectivityThresholds(2,2);
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

%% Prediction error and significance tests
for s = 1:length(data)
    % direction
    % only use units with minimum DSI and maximum OSI
    valid = find(data(s).DSI >= selectivityThresholds(1,1) & ...
        data(s).OSI <= selectivityThresholds(1,2) & ...
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
        data(s).DSI <= selectivityThresholds(2,2) & ...
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
    errDir_original = abs(dp - transVectors);
    ind = errDir_original > 180;
    errDir_original(ind) = 360 - errDir_original(ind);
    [errDir_original, optTransVects] = min(errDir_original, [], 2);
    figure
    histogram(optTransVects, .5:4.5)
    set(gca, 'XTick', 1:4, 'XTickLabel', ...
        {'advance','retreat','rise','fall'}, 'Box', 'off')
    title(sprintf('%s: Direction vectors', sets{s}))
    io.saveFigure(gcf, fPlots, sprintf('optimalVectors_%s_direction%s', ...
        sets{s}, suffix))
    % 2. orientation
    errOri_original = abs(op - [longVectors, latVectors]);
    ind = errOri_original > 90;
    errOri_original(ind) = 180 - errOri_original(ind);
    [errOri_original, optOriVects] = min(errOri_original, [], 2);
    figure
    histogram(optOriVects, .5:4.5)
    set(gca, 'XTick', 1:4, 'XTickLabel', ...
        {'H_{long}','V_{long}','H_{lat}','V_{lat}'}, 'Box', 'off')
    title(sprintf('%s: Orientation vectors', sets{s}))
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
    fprintf('      Mean: %.2f [%.2f-%.2f, %.2f-%.2f]\n', ...
        mean(errDir_original), ...
        prctile(mean(errDir_perm,1),[2.5 97.5]), ...
        prctile(mean(errDir_uni,1),[2.5 97.5]))
    fprintf('      Median: %.2f [%.2f-%.2f, %.2f-%.2f]\n', ...
        median(errDir_original), ...
        prctile(median(errDir_perm,1),[2.5 97.5]), ...
        prctile(median(errDir_uni,1),[2.5 97.5]))
    fprintf('    Orientation\n')
    fprintf('      Mean: %.2f [%.2f-%.2f, %.2f-%.2f]\n', ...
        mean(errOri_original), ...
        prctile(mean(errOri_perm,1),[2.5 97.5]), ...
        prctile(mean(errOri_uni,1),[2.5 97.5]))
    fprintf('      Median: %.2f [%.2f-%.2f, %.2f-%.2f]\n', ...
        median(errOri_original), ...
        prctile(median(errOri_perm,1),[2.5 97.5]), ...
        prctile(median(errOri_uni,1),[2.5 97.5]))
end
