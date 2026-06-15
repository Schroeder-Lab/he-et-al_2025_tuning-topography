function main_fitReceptiveFields_wEyePos(folders)

%% Parameters
sets = {'boutons', 'neurons'};

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

% for finding eye pos to gaze transformation
maxIter = 10;
minDelta = 1e-5;

%% Fit RFs and get explained variance
opts = struct();
% constraints on A = [a b; c d]:
% a <= 0 (if eye moves temporally (positive eye pos.), RF moves peripherally
% d >= 0 (if eye moves up, RF moves up as well)
opts.ABounds = [-0.5 -0.03; -0.5 0.5; -0.5 0.5; 0.03 0.5];
% horizontal/vertical eye movements has higher impact on horizontal/
% vertical RF movement:
% a <= |b| => a-b <= 0 and a+b <= 0
% d >= |c| => c-d <= 0 and -c-d <= 0
% -a >= d  => a+d <= 0
opts.linearA = [1 -1 0 0; 1 1 0 0; 0 0 1 -1; 0 0 -1 -1; 1 0 0 1];
opts.linearb = [0; 0; 0; 0; 0];

for s = 1:2 % boutons and neurons
    subjDirs = dir(fullfile(folders.data, sets{s}, 'SS*'));
    for subj = 1:length(subjDirs) % animals
        name = subjDirs(subj).name;
        dateDirs = dir(fullfile(folders.data, sets{s}, name, '2*'));
        for dt = 1:length(dateDirs) %dates
            date = dateDirs(dt).name;
            fprintf('\n\n%s  %s\n', name, date)
            f = fullfile(folders.data, sets{s}, name, date);
            % ignore session if visual noise stimulus was not presented or
            % pupil data is missing
            if ~isfile(fullfile(f, '_ss_sparseNoise.times.npy')) || ...
                    ~isfile(fullfile(f, 'eye.xyPos.npy'))
                continue
            end

            %% Load data
            caData = io.getCalciumData(f);
            stimData = io.getVisNoiseInfo(f);
            pupilData = io.getPupilData(f);

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

            % reshape stimulus frames to [time x px]; this represents a single
            % "stimulus block", i.e. the pixels to estimate a single time point of the
            % receptive field
            stim = reshape(stimMatrix, size(stimMatrix,1), []);

            % fill time gaps in stimulus with zeros
            [t_stim, stim] = stimuli.fillGapsInNoiseStim(t_stim, stim);

            % generate topelitz matrix, including time lags in RFs
            toeplitz = rf.makeStimToeplitz(stim, rfBins);

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

            % resample traces to stimulus frame times
            tBin_ca = median(diff(t_ca));
            numBins = round(tBin_stim / tBin_ca);
            zTraces = smoothdata(caTraces, 1, 'movmean', numBins, 'omitnan');
            zTraces = interp1(t_ca, zTraces, t_stim);

            % z-score neural response
            zTraces = (zTraces - mean(zTraces,1,'omitnan')) ./ ...
                std(zTraces,0,1,'omitnan');

            %% Prepare pupil data
            % estimate pupil position for each stimulus frame
            t_ind = pupilData.time > t_stim(1) - 1 & ...
                pupilData.time < t_stim(end) + 1;
            t_pupil = pupilData.time(t_ind);
            eyePos_cont = pupilData.xyPos(t_ind,:);

            eyePos = traces.getAlignedTraces(eyePos_cont, t_pupil, ...
                t_stim, [0 tBin_stim]); % [bins, frames, 2]
            eyePos = squeeze(mean(eyePos, 1, "omitnan")); % [frames, 2]

            % subtract median pupil position to get pupil shifts from
            % median
            defaultEyePos = median(eyePos, 1, "omitnan");
            eyePos = eyePos - defaultEyePos;

            % invert y-positions so that positive values reflect that pupil
            % is above median position
            eyePos(:,2) = -eyePos(:,2);

            %% Map RFs given current linear transform A (eye shift -> RF shift)
            A = [-0.04, 0; 0, 0.04];
            deltaShift = Inf;
            error = Inf;
            count = 0;

            while deltaShift > minDelta && count < maxIter

                fprintf('  Iteration %d\n\n', count + 1)

                % Determine shift of stimulus, which is in opposite
                % direction to shift of RFs
                stimShifts = eyePos * A'; % in vis. deg., [frames, 2]

                % create toeplitz matrix with shifted stimulus frames (consider
                % upscaling stimulus matrix)
                stim_shifted = stimuli.shiftNoiseFrames(stim, stimSize, ...
                    gridY, gridX, stimShifts);
                toeplitz_shifted = rf.makeStimToeplitz(stim_shifted, rfBins);

                % map pixel-RFs
                % rFields: [rows x columns x time x ON/OFF x units]
                rFields = rf.getReceptiveField( ...
                    zTraces, toeplitz_shifted, stimSize, rfBins, lambda);
                % NOTE: OFF subfield: positive pixel -> suppressed by black
                % negative pixel -> driven by black

                % fit Gaussians
                [rfGaussPars, fitGaussians, fitWeights, peakNoiseRatio, ...
                    bestSubFields, subFieldSigns] = ...
                    rf.fitAllRFs(rFields, rfBins, gridX, gridY, zTraces);

                % Optimize linear transform A (eye shift -> RF shift)
                % define a function that convolves current RF at a given
                % shift with the stimulus to produce the predicted neural
                % response
                opts.A0 = A;
                [AHat, out] = rf.fitFullEyeAffine_Astep(toeplitz, zTraces, ...
                    eyePos, rfGaussPars, bestSubFields, subFieldSigns, ...
                    gridX, gridY, rfBins, fitWeights, opts);

                % only update A if results better than previously
                if out.sse < error
                    error = out.sse;
                    deltaShift = mean((A - AHat).^2, "all");
                    A = AHat;
                else
                    break
                end
                count = count + 1;
                fprintf('  A = [%.2f, %.2f; %.2f, %.2f]\n\n', A')
            end

            % Using the best A (resulting in smallest error) and the RFs 
            % fitted on stimuli shifted according to A, predict responses 
            % and get EVs
            [predictions, EVs] = rf.makeAllNoisePredictions(zTraces, ...
                toeplitz_shifted, fitGaussians, fitWeights);

            %% Save results
            writeNPY(permute(rFields, [5 1 2 3 4]), fullfile(f, '_ss_rf_eye.maps.npy'))
            writeNPY(bestSubFields, fullfile(f, '_ss_rf_eye.bestSubField.npy'))
            writeNPY(subFieldSigns, fullfile(f, '_ss_rf_eye.subFieldSigns.npy'))
            writeNPY(rfGaussPars, fullfile(f, '_ss_rf_eye.gaussFitPars.npy'))
            writeNPY(peakNoiseRatio, fullfile(f, '_ss_rf_eye.peak.npy'))
            writeNPY(fitGaussians, fullfile(f, '_ss_rf_eye.gaussMask.npy'))
            writeNPY(fitWeights, fullfile(f, '_ss_rf_eye.gaussTimeWeights.npy'))
            writeNPY(EVs, fullfile(f, '_ss_rf_eye.explVars.npy'))
            writeNPY(predictions, fullfile(f, '_ss_rfPrediction_eye.traces.npy'))
            writeNPY(t_stim, fullfile(f, '_ss_rfPrediction_eye.timestamps.npy'))
            writeNPY(edges, fullfile(f, '_ss_rfDescr.edges.npy'))
            writeNPY(rfBins * tBin_stim, fullfile(f, '_ss_rfDescr.timestamps.npy'))

            writeNPY(eyePos, fullfile(f, '_ss_rfPrediction_eye.eyePos.npy'))
            writeNPY(A, fullfile(f, 'eyeToGaze.transformMatrix.npy'))
            writeNPY(defaultEyePos', fullfile(f, 'eyeToGaze.defaultEyePos.npy')) % in original pixel coordinates
        end
    end
end

%% Plot RFs
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
%             % ignore session if visual noise stimulus was not presented or
%             % pupil data is missing
%             if ~isfile(fullfile(f, '_ss_sparseNoise.times.npy')) || ...
%                     ~isfile(fullfile(f, 'eye.xyPos.npy'))
%                 continue
%             end
%             fPlots = fullfile(folders.plots, 'ReceptiveFields_eyePosCorrected', ...
%                 sets{s}, name, date);
%             if ~isfolder(fPlots)
%                 mkdir(fPlots)
%             end
% 
%             % load data
%             results = io.getNoiseRFFits(f, true);
% 
%             edges = results.edges;
%             gridW = diff(edges(1:2)) / size(results.maps,3);
%             gridH = -diff(edges(3:4)) / size(results.maps,2);
%             for iUnit = 1:size(results.gaussPars,1)
%                 % plot RF (if explained var. and peak-to-noise high enough
%                 if results.EV(iUnit) >= minEV_plot && ...
%                         results.peaks(iUnit) >= minPeak_plot
%                     % rf_tmp: [rows x columns x subfield]
%                     rf_tmp = squeeze(mean(results.maps(iUnit,:,:,:,:),4));
% 
%                     rf.plotRF(rf_tmp, results.gaussPars(iUnit,:), ...
%                         results.bestSubFields(iUnit), ...
%                         edges, edges, gridW, gridH)
% 
%                     sgtitle(sprintf('ROI %d (EV: %.3f, peak/noise: %.1f, %s)', ...
%                         iUnit, results.EV(iUnit), results.peaks(iUnit), ...
%                         RFtypes{results.bestSubFields(iUnit)}))
% 
%                     saveas(gcf, fullfile(fPlots, sprintf('Unit%03d.jpg', iUnit)));
%                     close gcf
%                 end
%             end
%         end
%     end
% end

%% Plot all RF outlines per dataset
for s = 1:2 % boutons and neurons
    fPlots = fullfile(folders.plots, 'ReceptiveFields_eyePosCorrected', sets{s});
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
            % ignore session if visual noise stimulus was not presented or
            % pupil data is missing
            if ~isfile(fullfile(f, '_ss_sparseNoise.times.npy')) || ...
                    ~isfile(fullfile(f, 'eye.xyPos.npy'))
                continue
            end

            data = io.getNoiseRFFits(f, true);

            rf.plotRFOutlines(data.gaussPars, data.EV, data.peaks, ...
                minEV, minPeak, [], data.edges, [])

            title(sprintf('%s %s', name, date))

            saveas(gcf, fullfile(fPlots, sprintf('RFs_%s_%s.jpg', name, date)));
            close gcf
        end
    end
end

%% Plot eye positions during RF mapping per dataset
for s = 1:2 % boutons and neurons
    fPlots = fullfile(folders.plots, 'ReceptiveFields_eyePosCorrected', sets{s});
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
            % ignore session if visual noise stimulus was not presented or
            % pupil data is missing
            if ~isfile(fullfile(f, '_ss_sparseNoise.times.npy')) || ...
                    ~isfile(fullfile(f, 'eye.xyPos.npy'))
                continue
            end

            data = io.getNoiseRFFits(f, true);
            A = data.A;

            gazePos = A * data.eyePos';

            figure
            subplot(2,1,1)
            scatter(gazePos(1,:), gazePos(2,:), 15, 'k', 'filled', 'o')
            axis image
            xlabel('\DeltaAzimuth (vis. deg.)')
            ylabel('\DeltaElevation (vis. deg.)')
            title(sprintf('%s %s', name, date))

            subplot(2,1,2)
            scatter(data.eyePos(:,1), data.eyePos(:,2), 15, 'k', 'filled', 'o')
            axis image
            set(gca, 'XDir', 'reverse')
            xlabel('\DeltaX (pixels)')
            ylabel('\DeltaY (pixels)')

            str = sprintf('A =\n[%.2f  %.2f\n %.2f  %.2f]', ...
                A(1,1), A(1,2), A(2,1), A(2,2));

            annotation('textbox', [0.9 0.3 0.1 0.1], ...
                'String', str, ...
                'LineStyle', 'none');

            saveas(gcf, fullfile(fPlots, sprintf('gazePos_%s_%s.jpg', name, date)));
            close gcf
        end
    end
end
