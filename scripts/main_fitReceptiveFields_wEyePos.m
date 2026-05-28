function main_fitReceptiveFields_wEyePos(folders)

%% Parameters
sets = {'boutons', 'neurons'};

% for receptive field fits
lambda = 0.002; % smoothing parameter for spatial RF
rf_timeLimits = [0.2 0.5]; % time range of stimulus before neural response 
                      % considered for RF
RFtypes = {'ON', 'OFF', 'ON+OFF'};

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
            stim = reshape(stimFrames, size(stimFrames,1), []);
            
            % fill time gaps in stimulus with zeros
            [t_stim, stim] = stimuli.fillGapsInNoiseStim(t_stim, stim);

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
            eyePos = pupilData.xyPos(t_ind,:);

            eyePos_frame = traces.getAlignedTraces(eyePos, t_pupil, ...
                t_stim, [0 tBin_stim]); % [bins, frames, 2]
            eyePos_frame = squeeze(mean(eyePos_frame, 1, "omitnan")); % [frames, 2]

            % subtract median pupil position to get pupil shifts from
            % median
            eyePos_frame = eyePos_frame - median(eyePos_frame, 1, "omitnan");
            
            %% Map RFs given current linear transform A (eye shift -> RF shift)
            A = zeros(2);

            stimShifts = eyePos_frame * A; % in vis. deg., [frames, 2]

            % create toeplitz matrix with shifted stimulus frames (consider
            % upscaling stimulus matrix)
            stim_shifted = stimuli.shiftNoiseFrames(stim, stimSize, ...
                gridY, gridX, stimShifts);
            toeplitz = rf.makeStimToeplitz(stim_shifted, t_stim, rfBins);

            % map pixel-RFs
            % rFields: [rows x columns x time x ON/OFF x units]
            rFields = rf.getReceptiveField( ...
                zTraces, toeplitz, stimSize, rfBins, lambda);
            % NOTE: OFF subfield: positive pixel -> suppressed by black
            % negative pixel -> driven by black

            % fit Gaussians
            [rfGaussPars, fitGaussians, fitWeights, peakNoiseRatio, ...
                bestSubFields, subFieldSigns] = ...
                rf.fitAllRFs(rFields, rfBins, gridX, gridY, zTraces);

            %% Optimize linear transform A (eye shift -> RF shift)
            % define a function that convolves current RF at a given
            % shift with the stimulus to produce the predicted neural
            % response

            % define a function that evaluates the prediction across all
            % neurons

            % optimize A by minimizing the error of all predicted neural
            % responses


        end
    end
end