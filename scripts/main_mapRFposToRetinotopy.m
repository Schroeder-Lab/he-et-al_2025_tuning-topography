% Find mapping between RF position and brain position (retinotopy). Then
% infer RF position of units where RF could not be mapped.

%% Folders
getFolders;

%% Parameters
minEV = 0.01;
minPeak = 3.5;

%% Add paths
addpath(genpath(fullfile(folders.tools, 'npy-matlab')))
addpath(fullfile(folders.repo))

%% Fit RFs and get cross-validated explained variance
sets = {'boutons', 'neurons'};
for s = 1:2 % boutons and neurons
    fPlots = fullfile(folders.plots, 'RetinotopicRFs', sets{s});
    if ~isfolder(fPlots)
        mkdir(fPlots)
    end

    subjDirs = dir(fullfile(folders.data, sets{s}, 'SS*'));
    for subj = 1:length(subjDirs) % animals
        name = subjDirs(subj).name;
        fprintf('%s\n', name)
        dateDirs = dir(fullfile(folders.data, sets{s}, name, '2*'));
        for dt = 1:length(dateDirs) %dates
            date = dateDirs(dt).name;
            fprintf('  %s\n', date)
            f = fullfile(folders.data, sets{s}, name, date);
            % ignore session if visual noise stimuolus was not present
            if ~isfile(fullfile(f, '_ss_sparseNoise.times.npy'))
                continue
            end

            % load data
            data = io.getRecordingInfo(f);
            brainPos = data.roiPositions; % (x,y,z) in microns
            brainPos(:,3) = [];
            data = io.getVisNoiseInfo(f);
            edges = data.edges;
            data = io.getRFFits;
            fits = data.fitParameters;
            rfPos = fits(:, [2 4]); % (azimuth, elevation) in visual degrees
            ev_rf = data.EV;

            % for all units with a significant RF, use linear regression to
            % model:
            % rf_x = a1 + b1 * brain_x + c1 * brain_y
            % rf_y = a2 + b2 * brain_x + c2 * brain_y
            valid = ev_rf >= minEV;
            fit_rf_x = fit(brainPos(valid,:), rfPos(valid,1), 'poly11', ...
                'Robust', 'Bisquare');
            fit_rf_y = fit(brainPos(valid,:), rfPos(valid,2), 'poly11', ...
                'Robust', 'Bisquare');

            % to directly check fitting, do those plots:
            figure
            plot(fit_rf_x, brainPos(valid,:), rfPos(valid,1))
            xlabel('Brain in X')
            ylabel('Brain in Y')
            zlabel('RF azimuth')
            title('Fit horizontal RF position')
            figure
            plot(fit_rf_y, brainPos(valid,:), rfPos(valid,2))
            xlabel('Brain in X')
            ylabel('Brain in Y')
            zlabel('RF elevation')
            title('Fit vertical RF position')

            % for all units, use linear model to predict RF position
            predict_rf_x = fit_rf_x(brainPos);
            predict_rf_y = fit_rf_y(brainPos);

            % save results
            writeNPY([predict_rf_x, predict_rf_y], ...
                fullfile(f, '_ss_rf.posRetinotopy.npy'))

            % plot
            figure('WindowState', 'maximized')
            tiledlayout(1,2)

            nexttile
            h = [0 0];
            plot([rfPos(valid,1)'; predict_rf_x(valid)'], ...
                [rfPos(valid,2)'; predict_rf_y(valid)'], 'k')
            hold on
            h(1) = plot(rfPos(valid,1), rfPos(valid,2), 'pentagram', ...
                'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'r', ...
                'MarkerSize', 10);
            h(2) = plot(predict_rf_x(valid), predict_rf_y(valid), 'o', ...
                'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'k');
            axis image
            set(gca, 'YDir', 'reverse', 'box', 'off')
            axis(edges)
            legend(h, 'Measured RF', 'Fitted RF')
            title('Measured vs fitted RF positions')

            nexttile
            hold on
            h = [0 0];
            h(1) = plot(predict_rf_x(~valid), predict_rf_y(~valid), 'o', ...
                'MarkerEdgeColor', 'none', 'MarkerFaceColor', [1 1 1].*0.75);
            h(1) = plot(predict_rf_x(valid), predict_rf_y(valid), 'o', ...
                'MarkerEdgeColor', 'none', 'MarkerFaceColor', [1 1 1].*0.25);
            axis image
            set(gca, 'YDir', 'reverse', 'box', 'off')
            legend(h, 'Unit without RF', 'Unit with RF')
            title('Fitted RF positions of all units')
            limits = [min(predict_rf_x) max(predict_rf_x) ...
                min(predict_rf_y) max(predict_rf_y)];
            limits([1 3]) = min([limits([1 3]); edges([1 3])], [], 1);
            limits([2 4]) = min([limits([2 4]); edges([2 4])], [], 1);
            axis limits
            if ~all(limits == edges)
                h(3) = plot(edges([1 1 2 2 1]), edges([3 4 4 3 3]), 'k:');
                legend(h, 'Unit without RF', 'Unit with RF', 'Monitor edge')
            end
            sgtitle(sprintf('%s %s', name, date), 'FontWeight', 'bold')
            
            saveas(gcf, fullfile(fPlots, sprintf('%s_%s.jpg', ...
                name, date)))
            close gcf
        end
    end
end