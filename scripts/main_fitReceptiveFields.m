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
% - best subfield is found by MSE of Gaussian fitted to RF map (either only
%   ON, only OFF, or both) (instead of comparing ratios of ON and OFF
%   peaks)
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
rf_timeLimits = [0.2 0.5]; % time range of stimulus before neural response 
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

titles = {'ON field','OFF field'};

sets = {'boutons', 'neurons'};

%% Fit RFs and get explained variance
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
            [stimMatrix, edges, gridX, gridY] = ...
                stimuli.getNoiseStimMatrix(stimData.edges, ...
                stimData.frames, stimData.stimOrder);
            stimSize = size(stimMatrix, [2 3]);
            t_stim = stimData.times;
            tBin_stim = median(diff(t_stim));
            rfBins = floor(rf_timeLimits(1) / tBin_stim) : ...
                ceil(rf_timeLimits(2) / tBin_stim) - 1;

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
            % reshape stimulus frames to [time x px]; this represents a single
            % "stimulus block", i.e. the pixels to estimate a single time point of the
            % receptive field
            stim = reshape(stimMatrix, size(stimMatrix,1), []);
            
            % fill time gaps in stimulus with zeros
            [t_stim, stim] = stimuli.fillGapsInNoiseStim(t_stim, stim);

            % generate toplitz matrix for stimulus
            toeplitz = rf.makeStimToeplitz(stim, rfBins);

            % resample neural response at stimulus times
            tBin_ca = median(diff(t_ca));
            numBins = round(tBin_stim / tBin_ca);
            zTraces = smoothdata(caTraces, 1, 'movmean', numBins, 'omitnan');
            zTraces = interp1(t_ca, zTraces, t_stim);
            % z-score neural response
            zTraces = (zTraces - mean(zTraces,1,'omitnan')) ./ ...
                std(zTraces,0,1,'omitnan');

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

            % fit Gaussian to mapped RFs
            [rfGaussPars, fitGaussians, fitWeights, peakNoiseRatio, ...
                bestSubFields, subFieldSigns] = ...
                rf.fitAllRFs(rFields, rfBins, gridX, gridY, zTraces);

            % predict responses and get EVs
            [predictions, EVs] = rf.makeAllNoisePredictions(zTraces, ...
                toeplitz, fitGaussians, fitWeights);

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
            writeNPY(t_stim, fullfile(f, '_ss_rfPrediction.timestamps.npy'))
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
            results = getNoiseRFFits(f);

            edges = results.edges;
            gridW = diff(edges(1:2)) / size(results.maps,3);
            gridH = -diff(edges(3:4)) / size(results.maps,2);
            for iUnit = 1:size(results.fitParameters,1)
                % plot RF (if explained var. and peak-to-noise high enough
                if results.EV(iUnit) >= minEV_plot && ...
                        results.peaks(iUnit) >= minPeak_plot
                    % rf_tmp: [rows x columns x subfield]
                    rf_tmp = squeeze(mean(results.maps(iUnit,:,:,:,:),4));

                    rf.plotRF(rf_tmp, results.fitParameters(iUnit,:), ...
                        results.bestSubFields(iUnit), ...
                        edges, edges, gridW, gridH)

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

            data = io.getNoiseRFFits(f);

            rf.plotRFOutlines(data.fitParameters, data.EV, data.peaks, ...
                minEV, minPeak, [], data.edges, [])

            title(sprintf('%s %s', name, date))

            saveas(gcf, fullfile(fPlots, sprintf('%s_%s.jpg', name, date)));
            close gcf
        end
    end
end
