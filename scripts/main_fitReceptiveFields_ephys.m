function main_fitReceptiveFields_ephys(folders)
% Given the spike responses to the visual noise stimulus, fit spatial
% receptive fields (ON and OFF fields).

%% Parameters
% for receptive field fits
lambda = 0.002; % smoothing parameter for spatial RF
rf_timeLimits = [0 0.2]; % time range of stimulus before neural response 
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

%% Fit RFs and get cross-validated explained variance
subjDirs = dir(fullfile(folders.data, 'ephys'));
subjDirs = subjDirs(~startsWith({subjDirs.name}, '.'));
for subj = 1:length(subjDirs) % animals
    name = subjDirs(subj).name;
    fprintf('%s\n', name)
    dateDirs = dir(fullfile(folders.data, 'ephys', name, '2*'));
    for dt = 1:length(dateDirs) %dates
        date = dateDirs(dt).name;
        fprintf('  %s\n', date)
        f = fullfile(folders.data, 'ephys', name, date);
        % ignore session if visual noise stimulus was not present
        if ~isfile(fullfile(f, '_ss_sparseNoise.times.npy'))
            continue
        end

        %% Load data
        spikeData = io.getEphysData(f);
        stimData = io.getVisNoiseInfo(f);
        % ignore spike data before/after visual noise stimulation
        t_ind = spikeData.times >= stimData.times(1) - 10 & ...
            spikeData.times <= stimData.times(end) + 10;
        spikeData.times = spikeData.times(t_ind);
        spikeData.clusters = spikeData.clusters(t_ind);

        %% Prepare stimulus data
        % edges: [left right top bottom] (above horizon: >0)
        [stimMatrix, edges, gridX, gridY] = ...
            stimuli.getNoiseStimMatrix(stimData.edges, ...
            stimData.frames, stimData.stimOrder);
        stimSize = size(stimMatrix, [2 3]);
        t_stim = stimData.times;
        tBin_stim = median(diff(t_stim));
        rfBins = floor(rf_timeLimits(1) / tBin_stim) : ...
            ceil(rf_timeLimits(2) / tBin_stim);

        %% Map RFs

        % generate toplitz matrix for stimulus
        [toeplitz, t_toeplitz] = ...
            rf.makeStimToeplitz(stimMatrix, t_stim, rfBins);

        % get firing rates aligned to stimulus times
        units = unique(spikeData.clusters);
        traces = nan(length(t_stim), length(units));
        for iUnit = 1:length(units)
            t = spikeData.times(spikeData.clusters == units(iUnit));
            [~, frameOfSpike] = events.alignData(t, ...
                t_stim, [0 tBin_stim]);
            traces(:,iUnit) = histcounts(frameOfSpike, ...
                (0:length(t_stim)) + 0.5);
        end

        % z-score neural response
        zTraces = (traces - mean(traces,1,'omitnan')) ./ ...
            std(traces,0,1,'omitnan');

        %--------------------------------------------------------------
        % Comment if RFs are already mapped and only Gaussian fit is
        % needed.
        % map RF
        % rFields: [rows x columns x time x ON/OFF x units]
        rFields = rf.getReceptiveField( ...
            zTraces, toeplitz, stimSize, rfBins, lambda);
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
        [rfGaussPars, fitGaussians, fitWeights, peakNoiseRatio, ...
            bestSubFields, subFieldSigns, predictions, EVs] = ...
            rf.fitAllRFs(rFields, rfBins, gridX, gridY, zTraces, toeplitz);

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
                    rf_tmp(:,:,2) = -rf_tmp(:,:,2);
                    mx = max(abs(rf_tmp),[],"all");
                    rfGaussPars = results.fitParameters(iUnit,:);
                    figure('Position', [75 195 1470 475])
                    for sf = 1:2
                        subplot(1,2,sf)
                        % STA
                        imagesc([edges(1)+gridW/2 edges(2)-gridW/2], ...
                            [edges(3)-gridH/2 edges(4)+gridH/2], ...
                            rf_tmp(:,:,sf),[-mx mx])
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
                        iUnit, results.EV(iUnit), results.peaks(iUnit), ...
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