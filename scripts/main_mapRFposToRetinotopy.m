function main_mapRFposToRetinotopy(folders)
% Find mapping between RF position and brain position (retinotopy). Then
% infer RF position of units where RF could not be mapped.

%% Parameters
minUnits = 10;
minEV = 0.01;
minPeak = 5;
% minEV = [0.01 0];
% minPeak = [6 12];
sets = {'boutons', 'neurons'};

%% Use linear regression to fit RF based on brain position 
% (one model for azimuth and elevation separately)
% for s = 1:2 % boutons and neurons
%     subjDirs = dir(fullfile(folders.data, sets{s}, 'SS*'));
%     for subj = 1:length(subjDirs) % animals
%         name = subjDirs(subj).name;
%         fprintf('%s\n', name)
%         dateDirs = dir(fullfile(folders.data, sets{s}, name, '2*'));
%         for dt = 1:length(dateDirs) %dates
%             date = dateDirs(dt).name;
%             fprintf('  %s\n', date)
%             f = fullfile(folders.data, sets{s}, name, date);
%             % ignore session if visual noise stimuolus was not present
%             if ~isfile(fullfile(f, '_ss_sparseNoise.times.npy'))
%                 continue
%             end
% 
%             % load data
%             data = io.getRecordingInfo(f);
%             brainPos = data.roiPositions; % (x,y,z) in microns
%             brainPos(:,3) = [];
%             data = io.getRFFits(f);
%             fits = data.fitParameters;
%             rfPos = fits(:, [2 4]); % (azimuth, elevation) in visual degrees
%             ev_rf = data.EV;
%             rf_peaks = data.peaks;
% 
%             % only consider significant RFs
%             valid = ev_rf >= minEV & rf_peaks >= minPeak;
%             % valid = (ev_rf > minEV(1) & rf_peaks > minPeak(1)) | ...
%             %     (ev_rf > minEV(2) & rf_peaks > minPeak(2));
%             if sum(valid) < minUnits
%                 continue
%             end
% 
%             % find outliers
%             outliers = false(length(valid),1);
%             outliers(valid) = any(isoutlier(rfPos(valid,:), ...
%                 "ThresholdFactor", 5),2);
%             valid = valid & ~outliers;
% 
%             ft = fittype( @(a1, b1, c1, a2, k2, x) ...
%                 rf.brainToRFPos(x, a1, b1, c1, a2, k2));
%             % ft = fittype('rf.brainToRFPos(brainPos, a1, b1, c1, a2, k2)', ...
%             %     'independent', 'brainPos');
%             fit_rf = fit(reshape(brainPos(valid,:),[],1), ...
%                 reshape(rfPos(valid,:),[],1), ft, ...
%                 'StartPoint', [0 1 1 0 1], ...
%                 'Weights', repmat(ev_rf(valid),2,1));
% 
%             % % for all units with a significant RF, use linear regression to
%             % % model:
%             % % rf_x = a1 + b1 * brain_x + c1 * brain_y
%             % % rf_y = a2 + b2 * brain_x + c2 * brain_y
%             % fit_rf_x = fit(brainPos(valid,:), rfPos(valid,1), 'poly11', ...
%             %     'Weights', ev_rf(valid));
%             % fit_rf_y = fit(brainPos(valid,:), rfPos(valid,2), 'poly11', ...
%             %     'Weights', ev_rf(valid));
% 
%             % % to directly check fitting, do those plots:
%             % figure
%             % plot(fit_rf_x, brainPos(valid,:), rfPos(valid,1))
%             % xlabel('Brain in X')
%             % ylabel('Brain in Y')
%             % zlabel('RF azimuth')
%             % title('Fit horizontal RF position')
%             % figure
%             % plot(fit_rf_y, brainPos(valid,:), rfPos(valid,2))
%             % xlabel('Brain in X')
%             % ylabel('Brain in Y')
%             % zlabel('RF elevation')
%             % title('Fit vertical RF position')
% 
%             predict_rf = fit_rf(reshape(brainPos,[],1));
%             predict_rf = reshape(predict_rf, [], 2);
% 
%             % predict_rf_x = fit_rf_x(brainPos);
%             % predict_rf_y = fit_rf_y(brainPos);
% 
%             % save results
%             writeNPY(outliers, fullfile(f, '_ss_rf.outliers.npy'))
%             writeNPY(predict_rf, fullfile(f, '_ss_rf.posRetinotopy.npy'))
%             % writeNPY([predict_rf_x, predict_rf_y], ...
%             %     fullfile(f, '_ss_rf.posRetinotopy.npy'))
%             writeNPY(coeffvalues(fit_rf), ...
%                 fullfile(f, "_ss_rfRetinotopy.model.npy"))
%             % writeNPY([coeffvalues(fit_rf_x); coeffvalues(fit_rf_y)], ...
%             %     fullfile(f, "_ss_rfRetinotopy.model.npy"))
%         end
%     end
% end

%% Plot retinotopic maps
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
            if ~isfile(fullfile(f, '_ss_rf.posRetinotopy.npy'))
                continue
            end

            % load data
            data = io.getRecordingInfo(f);
            brainPos = data.roiPositions; % (x,y,z) in microns
            brainPos(:,3) = [];
            data = io.getCalciumData(f);
            planes = data.planes;
            ids = data.ids;
            data = io.getVisNoiseInfo(f);
            edges = data.edges;
            data = io.getRFFits(f);
            fitPars = data.fitParameters;
            rfPos = fitPars(:, [2 4]); % (azimuth, elevation) in visual degrees
            ev_rf = data.EV;
            rf_peaks = data.peaks;
            outliers = readNPY(fullfile(f, '_ss_rf.outliers.npy'));
            pred = readNPY(fullfile(f, '_ss_rf.posRetinotopy.npy'));
            predict_rf_x = pred(:,1);
            predict_rf_y = pred(:,2);
            model = readNPY(fullfile(f, "_ss_rfRetinotopy.model.npy"));
            fit_rfPos = @(x) rf.brainToRFPos(x, model(1), model(2), ...
                model(3), model(4), model(5));
            % fit_x = @(x,y) model(1,1) + model(1,2).*x + model(1,3).*y;
            % fit_y = @(x,y) model(2,1) + model(2,2).*x + model(2,3).*y;
            
            % only consider significant RFs
            valid = ev_rf >= minEV & rf_peaks >= minPeak;
            valid = valid & ~outliers;
            % valid = (ev_rf > minEV(1) & rf_peaks > minPeak(1)) | ...
            %     (ev_rf > minEV(2) & rf_peaks > minPeak(2));
            if sum(valid) < minUnits
                continue
            end
            
            % for all units with a significant RF, plot RF position and
            % colour code (1) horizontal and (2) vertical brain position
            brainLimits = reshape([floor(min(brainPos(valid,:),[],1)); ...
                ceil(max(brainPos(valid,:),[],1))], [], 1)';
            brainRange = [diff(brainLimits(1:2)) diff(brainLimits(3:4))];
            visualLimits = reshape([floor(min(rfPos(valid,:),[],1)); ...
                ceil(max(rfPos(valid,:),[],1))], [], 1)';
            retinotopyRange = [floor(min(predict_rf_x)), ceil(max(predict_rf_x)), ...
                floor(min(predict_rf_y)), ceil(max(predict_rf_y))];
            retinotopyRange([1 3]) = min(retinotopyRange([1 3]), ...
                visualLimits([1 3]));
            retinotopyRange([2 4]) = max(retinotopyRange([2 4]), ...
                visualLimits([2 4]));
            [x, y] = meshgrid(linspace(brainLimits(1)-brainRange(1), ...
                brainLimits(2)+brainRange(1), 60), ...
                linspace(brainLimits(3)-brainRange(2), ...
                brainLimits(4)+brainRange(2), 60));
            rfXY = fit_rfPos([x(:); y(:)]);
            rfX = reshape(rfXY(1:3600),60,60);
            rfY = reshape(rfXY(3601:end),60,60);
            % rfX = fit_x(x,y);
            % rfY = fit_y(x,y);
            plotLimits = brainLimits + 0.1.*[-brainRange(1), brainRange(1), ...
                -brainRange(2), brainRange(2)];
            contourLimits = fit_rfPos([plotLimits([1 1 2 2]), ...
                plotLimits([3 4 3 4])]);

            figure('Position', [80 255 1105 420])
            tiledlayout(1,2)

            nexttile
            hold on
            % contourLimits = fit_x(plotLimits([1 1 2 2]), plotLimits([3 4 3 4]));
            cLimits = [min([contourLimits(1:4) visualLimits(1)]), ...
                max([contourLimits(1:4) visualLimits(2)])];
            [~,c] = contourf(x, y, rfX, floor(cLimits(1)):floor(cLimits(2)));
            clim(cLimits)
            c.LineStyle = "none";
            scatter(brainPos(valid,1), brainPos(valid,2), [], ...
                rfPos(valid,1), "filled", 'MarkerEdgeColor', 'k')
            colormap turbo
            c = colorbar;
            c.Label.String = 'RF azimuth';
            axis image
            axis(plotLimits)
            set(gca, "Box", "off", "YDir", "reverse")
            xlabel('Brain ML (\mum)')
            ylabel('Brain AP (\mum)')
            title('Brain position vs RF azimuth')

            nexttile
            hold on
            % contourLimits = fit_y(plotLimits([1 1 2 2]), plotLimits([3 4 3 4]));
            cLimits = [min([contourLimits(5:end) visualLimits(3)]), ...
                max([contourLimits(5:end) visualLimits(4)])];
            [~,c] = contourf(x, y, rfY, floor(cLimits(1)):floor(cLimits(2)));
            clim(cLimits)
            c.LineStyle = "none";
            scatter(brainPos(valid,1), brainPos(valid,2), [], ...
                rfPos(valid,2), "filled", 'MarkerEdgeColor', 'k')
            colormap turbo
            c = colorbar;
            c.Label.String = 'RF elevation';
            axis image
            axis(plotLimits)
            set(gca, "Box", "off", "YDir", "reverse")
            xlabel('Brain ML (\mum)')
            ylabel('Brain AP (\mum)')
            title('Brain position vs RF elevation')
            sgtitle(sprintf('%s %s', name, date), 'FontWeight', 'bold')
            saveas(gcf, fullfile(fPlots, sprintf('%s_%s_RF-brain-position.jpg', ...
                name, date)))
            close gcf

            % plot fits
            figure('WindowState', 'maximized')
            drawnow
            sz = get(gcf, 'Position');
            set(gcf, 'Position', [sz(1) 2*sz(2) sz(3) round((sz(4)-sz(2))/2)])
            tiledlayout(1,3)

            nexttile
            hold on
            plot(rfPos(valid,1), rfPos(valid,2), 'or', 'LineWidth', 1.5)
            plot(rfPos(outliers,1), rfPos(outliers,2), 'xb', ...
                'LineWidth', 1.5)
            axis image
            ax = axis;
            ax([1 3]) = min(ax([1 3]), retinotopyRange([1 3]));
            ax([2 4]) = max(ax([2 4]), retinotopyRange([2 4]));
            axis(ax)
            set(gca, 'box', 'off')
            legend('Valid RFs', 'Outliers')
            xlabel('Azimuth (deg)')
            ylabel('Elevation (deg)')
            title('Outlier RF positions')

            nexttile
            h = {};
            plot([rfPos(valid,1)'; predict_rf_x(valid)'], ...
                [rfPos(valid,2)'; predict_rf_y(valid)'], 'k')
            hold on
            h{1} = plot(rfPos(valid,1), rfPos(valid,2), 'or', ...
                'LineWidth', 1.5);
            h{2} = plot(predict_rf_x(valid), predict_rf_y(valid), 'o', ...
                'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'k');
            dtRows = [dataTipTextRow("Plane", planes(valid)),...
                dataTipTextRow("ID", ids(valid)), ...
                dataTipTextRow("RF EV", ev_rf(valid)), ...
                dataTipTextRow("RF peak", rf_peaks(valid))];
            h{1}.DataTipTemplate.DataTipRows(3:6) = dtRows;
            axis image
            axis(retinotopyRange)
            set(gca, 'box', 'off')
            legend(cat(2,h{:}), 'Measured RF', 'Fitted RF')
            xlabel('Azimuth (deg)')
            ylabel('Elevation (deg)')
            title('Measured vs fitted RF positions')

            nexttile
            hold on
            h = [0 0 0];
            h(1) = plot(predict_rf_x(~valid), predict_rf_y(~valid), 'o', ...
                'MarkerEdgeColor', 'none', 'MarkerFaceColor', [1 1 1].*0.75);
            h(2) = plot(predict_rf_x(valid), predict_rf_y(valid), 'o', ...
                'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'k');
            h(3) = plot(edges([1 1 2 2 1]), edges([3 4 4 3 3]), 'k:');
            axis image
            axis(retinotopyRange)
            set(gca, 'box', 'off')
            legend(h, 'Unit without RF', 'Unit with RF', 'Monitor edge')
            xlabel('Azimuth (deg)')
            ylabel('Elevation (deg)')
            title('Fitted RF positions of all units')
            sgtitle(sprintf('%s %s', name, date), 'FontWeight', 'bold')
            
            saveas(gcf, fullfile(fPlots, sprintf('%s_%s_RF_mapped-fitted.jpg', ...
                name, date)))
            saveas(gcf, fullfile(fPlots, sprintf('%s_%s_RF_mapped-fitted', ...
                name, date)))
            close gcf
        end
    end
end