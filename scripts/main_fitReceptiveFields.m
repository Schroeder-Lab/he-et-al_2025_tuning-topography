function main_fitReceptiveFields(folders)
% Given the calcium responses to the visual noise stimulus, fit spatial
% receptive fields (ON and OFF fields).

% Differences to method used for saccade paper:
% - Removing initial decay in calcium traces: fit double-exponential to
%   whole trace (not just start)
% - RF (linear fit) is only smoothed in space, not in time (for saccade
%   paper also in time)
% - no crossvalidation
% - lambda is fixed
% - spatiotemporal RF is not interpolated before fitting the 2D Gaussian
% - explained variance determined from prediction based on Gaussian fit
%   (not from cross-validation)

%% Parameters
% for correcting baseline drifts of calcium traces at start of experiments
win_decay = 20; % in s, window to test whether baseline is higher than normal
thresh_decay = 1.5; % in std, threshold for decay

% to high-pass filter traces
smoothWin = 20; % in s

% for receptive field fits
lambda = 0.002; % smoothing parameter for spatial RF
rf_timeLimits = [0.2 0.4]; % time range of stimulus before neural response 
                      % considered for RF
RFtypes = {'ON', 'OFF', 'ON+OFF'};

% for evaluation of receptive fields (significance/goodness)
minEV = 0.01; % minimum explained variance to plot RF
minPeak = 5; % minimum peak of RF (compared to noise) to plot RF
minEV_plot = 0.005; % minimum explained variance to plot RF
minPeak_plot = 3.5; % minimum peak of RF (compared to noise) to plot RF

% plotting RFs
% colormaps
[cm_ON, cm_OFF] = colmaps.getRFMaps;
cms = cat(3, cm_ON, cm_OFF);
ellipse_x = linspace(-pi, pi, 100);

titles = {'ON field','OFF field'};

sets = {'boutons', 'neurons'};

%% Fit RFs and get cross-validated explained variance
for s = 1:2 % boutons and neurons
    subjDirs = dir(fullfile(folders.data, sets{s}, 'SS*'));
    for subj = 1:length(subjDirs) % animals
        name = subjDirs(subj).name;
        fprintf('%s\n', name)
        dateDirs = dir(fullfile(folders.data, sets{s}, name, '2*'));
        for dt = 1:length(dateDirs) %dates
            date = dateDirs(dt).name;
            fprintf('  %s\n', date)
            f = fullfile(folders.data, sets{s}, name, date);
            % ignore session if visual noise stimulus was not present
            if ~isfile(fullfile(f, '_ss_sparseNoise.times.npy'))
                continue
            end

            %% Load data
            caData = io.getCalciumData(f);
            stimData = io.getVisNoiseInfo(f);

            %% Prepare stimulus data
            % edges: [left right top bottom] (above horizon: >0)
            edges = double([stimData.edges([1 2]), -stimData.edges([3 4])]);
            gridW = diff(edges(1:2)) / size(stimData.frames,3);
            gridH = -diff(edges(3:4)) / size(stimData.frames,2);
            % ignore pixels in ipsilateral (right) hemifield
            if edges(1) * edges(2) < 0
                % determine right edge of all pixel columns
                rightEdges = edges(1) + ...
                    (1:size(stimData.frames,3)) .* gridW;
                validPix = find(rightEdges <= 0);
                stimData.frames = stimData.frames(:,:,validPix);
                edges(2) = rightEdges(validPix(end));
            end
            stimMatrix = stimData.frames(stimData.stimOrder,:,:);
            gridX = linspace(edges(1)+0.5*gridW, edges(2)-0.5*gridW, ...
                size(stimMatrix,3));
            gridY = linspace(edges(3)-0.5*gridH, edges(4)+0.5*gridH, ...
                size(stimMatrix,2));
            stimSize = size(stimMatrix, [2 3]);
            t_stim = stimData.times;
            tBin_stim = median(diff(t_stim));
            rfBins = floor(rf_timeLimits(1) / tBin_stim) : ...
                ceil(rf_timeLimits(2) / tBin_stim);

            %% Prepare calcium traces
            % ignore data before/after visual noise stimulation
            t_ind = caData.time > t_stim(1) - 10 & ...
                caData.time < t_stim(end) + 10;
            caTraces = caData.traces(t_ind,:);
            t_ca = caData.time(t_ind);
            % interpolate calcium traces to align all to same time
            [caTraces, t_ca] = traces.alignSampling(caTraces, t_ca, ...
                caData.planes, caData.delays);

            % remove strong baseline decay at start of experiment
            caTraces = traces.removeInitialDecay(caTraces, t_ca, ...
                win_decay, thresh_decay);

            % high-pass filter traces: remove smoothed traces
            caTraces = traces.highPassFilter(caTraces, t_ca, smoothWin);

            %% Map RFs

            % generate toplitz matrix for stimulus
            [toeplitz, t_toeplitz] = ...
                rf.makeStimToeplitz(stimMatrix, t_stim, rfBins);
            %--------------------------------------------------------------
            % Comment if RFs are already mapped and only Gaussian fit is
            % needed.
            % map RF
            % rFields: [rows x columns x time x ON/OFF x units]
            rFields = rf.getReceptiveField( ...
                caTraces, t_ca, toeplitz, t_toeplitz, stimSize, ...
                rfBins, lambda);
            % NOTE: OFF subfield: positive pixel -> suppressed by black
            % negative pixel -> driven by black

            %--------------------------------------------------------------
            % % Uncomment if RFs have already been mapped, and only Gaussian
            % % fit should be performed.
            % results = io.getRFFits(f);
            % rFields = results.maps;
            % rFields = permute(rFields, [2:5 1]);
            %--------------------------------------------------------------

            % fit Gaussian
            % parameters of fitted Gaussian:
            % [amplitude, xCenter, xStd, yCenter, yStd, rotation]
            rfGaussPars = NaN(size(caTraces,2), 6); 
            % fitted 2D Gaussian map (one for ON and OFF)
            fitGaussians = NaN(size(caTraces,2), size(rFields,1), size(rFields,2), 2);
            fitWeights = NaN(size(caTraces,2), length(rfBins));
            peakNoiseRatio = NaN(size(caTraces,2), 1);
            bestSubFields = NaN(size(caTraces,2), 1);
            subFieldSigns = NaN(size(caTraces,2), 2);
            predictions = NaN(length(t_toeplitz), size(caTraces,2));
            EVs = NaN(size(caTraces,2), 1);
            for iUnit = 1:size(rfGaussPars,1)
                % rfield: [rows x cols x t x ON/OFF]
                rfield = rFields(:,:,:,:,iUnit);
                % invert polarity of OFF field so that positive values
                % mean: unit is driven by black square
                rfield(:,:,:,2) = -rfield(:,:,:,2);
                % find best subfield (combination): find whether Gaussian 
                % is best fit to only ON, only OFF, or ON-OFF subfields, 
                % and what optimal sign of each subfield is
                % average across time
                rf_tmp = squeeze(mean(rfield,3));

                [fitRFs, RFsigns, MSEs] = rf.findRFGaussianMask(rf_tmp);
                [~, bestSubField] = min(MSEs);
                bestSubFields(iUnit) = bestSubField;
                subFieldSigns(iUnit,:) = RFsigns;
                fitGaussians(iUnit,:,:,:) = fitRFs(:,:,:,bestSubField);

                if bestSubField < 3
                    rf_sub = rf_tmp(:,:,bestSubField) .* RFsigns(bestSubField);
                else
                    rf_sub = (rf_tmp(:,:,1) .* RFsigns(1) + ...
                        rf_tmp(:,:,2) .* RFsigns(2)) ./ 2;
                end

                % fit Gaussian
                [rfGaussPars(iUnit,:), rf_gauss] = rf.fit2dGaussRF(rf_sub, false, ...
                    gridX, gridY);
                % mirror RF orientation to account for flipped y-axis direction
                % (top is positive)
                rfGaussPars(iUnit,6) = -rfGaussPars(iUnit,6);

                % subtract Gaussian from original RF map
                noise = rf_sub - rf_gauss;
                % distance of peak from noise
                peakNoiseRatio(iUnit) = rfGaussPars(iUnit,1) / std(noise(:));

                % predict response from RF
                % amplitudes (weights) of spatial RF across time span of RF
                weights = reshape(fitRFs(:,:,:,bestSubField), [], 1) \ ...
                    reshape(permute(rfield, [1 2 4 3]), [], size(rfield,3));
                % generate spatio-temporal RF from fitted Gaussian and
                % temporal weights
                spatTempMask = reshape(fitRFs(:,:,:,bestSubField), [], 1) * ...
                    weights; % [pix x t]
                spatTempMask = permute(reshape(spatTempMask, size(fitRFs,1), ...
                    size(fitRFs,2), 2, length(weights)), [1 2 4 3]); % [rows x cols x t x ON/OFF]
                % predict calcium trace based on generated spatio-temporal
                % RF
                [predictions(:, iUnit), EVs(iUnit)] = ...
                    rf.predictFromRF(caTraces(:,iUnit), t_ca, toeplitz, ...
                    t_toeplitz, spatTempMask);
                fitWeights(iUnit,:) = weights;
            end

            % save results
            writeNPY(permute(rFields, [5 1 2 3 4]), fullfile(f, '_ss_rf.maps.npy'))
            writeNPY(bestSubFields, fullfile(f, '_ss_rf.bestSubField.npy'))
            writeNPY(subFieldSigns, fullfile(f, '_ss_rf.subFieldSigns.npy'))
            writeNPY(rfGaussPars, fullfile(f, '_ss_rf.gaussFitPars.npy'))
            writeNPY(peakNoiseRatio, fullfile(f, '_ss_rf.peak.npy'))
            writeNPY(fitGaussians, fullfile(f, '_ss_rf.gaussMask.npy'))
            writeNPY(fitWeights, fullfile(f, '_ss_rf.gaussTimeWeights.npy'))
            writeNPY(EVs, fullfile(f, '_ss_rf.explVars.npy'))
            writeNPY(predictions, fullfile(f, '_ss_rfPrediction.traces.npy'))
            writeNPY(t_toeplitz, fullfile(f, '_ss_rfPrediction.timestamps.npy'))
            writeNPY(edges, fullfile(f, '_ss_rfDescr.edges.npy'))
            writeNPY(rfBins * tBin_stim, fullfile(f, '_ss_rfDescr.timestamps.npy'))
        end
    end
end

%% Plot RFs
for s = 1:2 % boutons and neurons
    subjDirs = dir(fullfile(folders.data, sets{s}, 'SS*'));
    for subj = 1:length(subjDirs) % animals
        name = subjDirs(subj).name;
        fprintf('%s\n', name)
        dateDirs = dir(fullfile(folders.data, sets{s}, name, '2*'));
        for dt = 1:length(dateDirs) %dates
            date = dateDirs(dt).name;
            fprintf('  %s\n', date)
            f = fullfile(folders.data, sets{s}, name, date);
            % ignore session if visual noise stimulus was not present
            if ~isfile(fullfile(f, '_ss_sparseNoise.times.npy'))
                continue
            end
            fPlots = fullfile(folders.plots, 'ReceptiveFields', ...
                sets{s}, name, date);
            if ~isfolder(fPlots)
                mkdir(fPlots)
            end

            % load data
            results = io.getRFFits(f);
            edges = results.edges;
            gridW = diff(edges(1:2)) / size(results.maps,3);
            gridH = -diff(edges(3:4)) / size(results.maps,2);
            for iUnit = 1:size(results.fitParameters,1)
                % plot RF (if explained var. and peak-to-noise high enough
                if results.EV(iUnit) >= minEV_plot && ...
                        results.peaks(iUnit) >= minPeak_plot
                    line = ':';
                    if results.EV(iUnit) >= minEV && ...
                            results.peaks(iUnit) >= minPeak
                        line = '-';
                    end
                    % rf_tmp: [rows x columns x subfield]
                    rf_tmp = squeeze(mean(results.maps(iUnit,:,:,:,:),4));
                    rfGaussPars = results.fitParameters(iUnit,:);
                    figure('Position', [75 195 1470 475])
                    for sf = 1:2
                        subplot(1,2,sf)
                        % STA
                        imagesc([edges(1)+gridW/2 edges(2)-gridW/2], ...
                            [edges(3)-gridH/2 edges(4)+gridH/2], ...
                            rf_tmp(:,:,sf) .* ...
                            results.subFieldSigns(iUnit,sf),[-mx mx])
                        hold on
                        if ismember(results.bestSubFields(iUnit), [sf 3])
                            % ellipse at 2 STD (x and y), not rotated, not shifted
                            x = rfGaussPars(3) * cos(ellipse_x);
                            y = rfGaussPars(5) * sin(ellipse_x);
                            % rotate and shift ellipse
                            x_rot = rfGaussPars(2) + ...
                                x .* cos(rfGaussPars(6)) - ...
                                y .* sin(rfGaussPars(6));
                            y_rot = rfGaussPars(4) + ...
                                x .* sin(rfGaussPars(6)) + ...
                                y .* cos(rfGaussPars(6));
                            n = x_rot < edges(1) | x_rot > edges(2) | ...
                                y_rot > edges(3) | y_rot < edges(4);
                            x_rot(n) = NaN;
                            y_rot(n) = NaN;
                            plot(x_rot, y_rot, ['k' line], 'LineWidth', 2)
                        end
                        axis image
                        set(gca, 'box', 'off', 'YDir', 'normal')
                        xlim(edges([1 2]))
                        ylim(edges([4 3])) 
                        colormap(gca, cms(:,:,sf))
                        title(titles{sf})
                        colorbar
                    end
                    sgtitle(sprintf('ROI %d (EV: %.3f, peak/noise: %.1f, %s)', ...
                        iUnit, results.EV(iUnit), results.peak(iUnit), ...
                        RFtypes{results.bestSubFields(iUnit)}))

                    saveas(gcf, fullfile(fPlots, sprintf('Unit%03d.jpg', iUnit)));
                    close gcf
                end
            end
        end
    end
end

%% Plot all RF outlines per dataset
for s = 1:2 % boutons and neurons
    fPlots = fullfile(folders.plots, 'ReceptiveFields', sets{s});
    if ~isfolder(fPlots)
        mkdir(fPlots)
    end

    subjDirs = dir(fullfile(folders.data, sets{s}, 'SS*'));
    for subj = 1:length(subjDirs) % animals
        name = subjDirs(subj).name;
        dateDirs = dir(fullfile(folders.data, sets{s}, name, '2*'));
        for dt = 1:length(dateDirs) %dates
            date = dateDirs(dt).name;
            f = fullfile(folders.data, sets{s}, name, date);
            % ignore session if visual noise stimulus was not present
            if ~isfile(fullfile(f, '_ss_sparseNoise.times.npy'))
                continue
            end

            caData = io.getRFFits(f);
            rfGaussPars = caData.fitParameters;
            EVs = caData.EV;
            peakNoiseRatio = caData.peaks;
            figure
            hold on
            for iUnit = 1:size(rfGaussPars,1)
                if EVs(iUnit) < minEV || peakNoiseRatio(iUnit) < minPeak
                    continue
                end
                % ellipse at 2 STD (x and y), not rotated, not shifted
                x = rfGaussPars(iUnit,3) * cos(ellipse_x);
                y = rfGaussPars(iUnit,5) * sin(ellipse_x);
                % rotate and shift ellipse
                x_rot = rfGaussPars(iUnit,2) + ...
                    x .* cos(rfGaussPars(iUnit,6)) - ...
                    y .* sin(rfGaussPars(iUnit,6));
                y_rot = rfGaussPars(iUnit,4) + ...
                    x .* sin(rfGaussPars(iUnit,6)) + ...
                    y .* cos(rfGaussPars(iUnit,6));
                plot(x_rot, y_rot, 'k')
            end
            axis image
            axis(caData.edges([1 2 4 3]))
            xlabel('Azimuth (visual degrees)')
            ylabel('Elevation (visual degrees)')
            title(sprintf('%s %s', name, date))

            saveas(gcf, fullfile(fPlots, sprintf('%s_%s.jpg', name, date)));
            close gcf
        end
    end
end

%% Plot example receptive fields for paper
% name = 'SS048';
% date = '2015-09-26';
% ROIs = [29 31 71 137];
% stimAzimuth = -110;
% stimElevation = 15;
% stimSigma = 9;
% RFstd = 1.5;
% RFheight = normpdf(RFstd) / normpdf(0);
% folder = fullfile(folderBase, name, date);
% 
% stimPos = readNPY(fullfile(folder, '_ss_sparseNoiseArea.edges.npy')); % [left right top bottom]
% rFields = readNPY(fullfile(folder, '_ss_rf.maps.npy'));
% fitPars = readNPY(fullfile(folder, '_ss_rf.gaussFitPars.npy'));
% subFields = readcell(fullfile(folder, '_ss_rf.type.csv'));
% 
% stimPos(3:4) = -stimPos(3:4);
% squW = diff(stimPos(1:2)) / size(rFields,3);
% squH = -diff(stimPos(3:4)) / size(rFields,2);
% [x0, y0] = meshgrid(linspace(stimPos(1), stimPos(2), 100), ...
%     flip(linspace(stimPos(4), stimPos(3), 100)));
% for ex = 1:length(ROIs)
%     rfield = squeeze(rFields(ROIs(ex),:,:,:,:)); % [rows x cols x t x ON/OFF]
%     rfield(:,:,:,2) = -rfield(:,:,:,2);
%     [mx,mxTime] = max(max(reshape(permute(abs(rfield), [1 2 4 3]), ...
%         [], size(rfield,3)), [], 1));
%     fitRF = whiteNoise.D2GaussFunctionRot(fitPars(ROIs(ex),:), cat(3, x0, y0));
%     figure('Position', [75 195 1470 475])
%     for f = 1:2
%         subplot(1,2,f)
%         imagesc([stimPos(1)+squW/2 stimPos(2)-squW/2], ...
%             [stimPos(3)-squH/2 stimPos(4)+squH/2], ...
%             rfield(:,:,mxTime,f),[-mx mx])
%         hold on
%         if any(strcmp(subFields{ROIs(ex)}, [RFtypes(f), RFtypes(3)]))
%             contour(x0, y0, fitRF, [1 1] .* (fitPars(ROIs(ex),1) * RFheight), ...
%                 'k', 'LineWidth', 2)
%         end
%         rectangle('Position', [stimAzimuth-stimSigma ...
%             stimElevation-stimSigma 2*stimSigma 2*stimSigma], ...
%             'Curvature', 1, 'LineStyle', '--', 'LineWidth', 2)
%         axis image
%         set(gca, 'box', 'off', 'XTick', round(stimPos(1:2)), ...
%             'YTick', [stimPos(4) 0 stimPos(3)], 'YDir', 'normal', ...
%             'YTickLabel', [stimPos(4) 0 stimPos(3)])
%         xlim(round([stimPos(1) stimPos(2)]))
%         colormap(gca, cms(:,:,f))
%         sgtitle(sprintf('ROI %d', ROIs(ex)))
%         colorbar
%     end
% 
%     filename = fullfile('C:\Users\Sylvia\OneDrive - University of Sussex\Paper\Plots\Fig 1', ...
%         sprintf('RF%03d.eps', ROIs(ex)));
%     myprint
%     close gcf
% end