function main_determineSCdepth(folders)

%% Parameters
smoothLFP = {5 21};
% stimWin = [0.02 .15];
stimWin = [-.05 .15];
baseWin = [-.05 0];
surfaceAmplitude = 0.25;

colors = colmaps.getBlueWhiteRedMap(255);
plotFolder = '07_SCDepthOnEphysProbe';
fPlots = fullfile(folders.results, plotFolder);
if ~isfolder(fPlots)
    mkdir(fPlots)
end

%% Determine SC depth from LFP
subjDirs = dir(fullfile(folders.data, 'ephys'));
subjDirs = subjDirs(~matches({subjDirs.name}, [".",".."]));
for subj = 1:length(subjDirs) % animals
    name = subjDirs(subj).name;
    dateDirs = dir(fullfile(folders.data, 'ephys', name, '2*'));
    for dt = 1:length(dateDirs) %dates
        date = dateDirs(dt).name;
        fprintf('%s %s\n', name, date)
        f = fullfile(folders.data, 'ephys', name, date);
        f_raw = fullfile(folders.dataRawEphys, name, date, 'exp01', ...
            sprintf('%s_%s_g0', name, date), ...
            sprintf('%s_%s_g0_imec0', name, date));
        
        %% Load data
        fprintf('  Load and prep LFP...\n')
        % load visual noise data
        stim = io.getVisNoiseInfo(f);
        if isempty(stim)
            fprintf('  NO VISUAL NOISE DATA. DISREGARD DATASET!\n')
            continue
        end

        % load LFP meta data
        meta = io.getLFPMetaData(fullfile(f_raw, ...
            sprintf('%s_%s_g0_t0.imec0.lf.meta', name, date)));
        channelLabels = readNPY(fullfile(f_raw, 'ibl_sorter', ...
            'channel_labels.npy'));
        channelCoords = readNPY(fullfile(f, 'channels.localCoordinates.npy'));
        channelDepths = channelCoords(1:2:end,2);
        chanDistance = diff(channelDepths(1:2)) / 1000; % millimeters

        % load LFP data
        fileLFP = dir(fullfile(f_raw, ...
            sprintf('%s_%s_g0_t0.imec0.lf.bin', name, date)));
        lfp = memmapfile(fullfile(f_raw, fileLFP.name), ...
            'Format',  {'int16', ...
            [meta.numChans fileLFP.bytes/(meta.numChans * 2)], 'x'}); % [channels x time]
        t_lfp = (1 : fileLFP.bytes/(meta.numChans * 2))' ./ meta.samplingRate;

        %% Prepare data
        % restrict LFP to stimulus period
        start = find(t_lfp > stim.times(1)-5, 1);
        stop = find(t_lfp > stim.times(end)+1, 1);
        % last channel is external input -> not LFP
        lfp = lfp.Data.x(1:end-1,start:stop); 
        t_lfp = t_lfp(start:stop);

        % conversion to Volt explained in
        % https://billkarsh.github.io/SpikeGLX/Support/Metadata_3A.html
        lfp = double(lfp) .* meta.Vmax  .* 1000 ./ meta.Imax ./ meta.lfpGain; % mV

        % subtract median across time from each channel
        lfp = lfp - median(lfp,2);

        % subtract median lfp trace (across channels) from each channel
        lfp_median = median(lfp, 1);
        lfp = lfp - lfp_median;

        % % Test whether filtering is necessary
        % [pxx, freq] = pwelch(lfp_mean, 512, 256, 1024, meta.samplingRate);
        % pxx_log = log10(pxx);
        % [~, ind50] = min(abs(freq - 50));
        % pow50 = pxx_log(ind50);
        % [~, ind40] = min(abs(freq - 40));
        % pow40 = pxx_log(ind40);
        % if pow50 - pow40 > 0.1
        %     f1 = figure;
        %     plot(t_lfp(150000:155000), lfp_mean(150000:155000))
        %     hold on
        %     f2 = figure;
        %     plot(freq, pxx)
        %     hold on
        % 
        %     % Remove 50Hz line noise
        %     fprintf('  Filter LFP...\n')
        %     notchFilter50 = designfilt('bandstopiir', ...
        %         'PassbandFrequency1',48, ...
        %         'StopbandFrequency1',49.5, ...
        %         'StopbandFrequency2',50.5, ...
        %         'PassbandFrequency2',52, ...
        %         'PassbandRipple1',1, ...
        %         'StopbandAttenuation',45, ...
        %         'PassbandRipple2',1, ...
        %         'SampleRate',meta.samplingRate);
        %     notchFilter100 = designfilt('bandstopiir', ...
        %         'PassbandFrequency1',95, ...
        %         'StopbandFrequency1',99, ...
        %         'StopbandFrequency2',101, ...
        %         'PassbandFrequency2',105, ...
        %         'PassbandRipple1',1, ...
        %         'StopbandAttenuation',30, ...
        %         'PassbandRipple2',1, ...
        %         'SampleRate',meta.samplingRate);
        % 
        %     % notchFilter6 = designfilt('bandstopiir', ...
        %     %     'PassbandFrequency1',5, ...
        %     %     'StopbandFrequency1',6, ...
        %     %     'StopbandFrequency2',7, ...
        %     %     'PassbandFrequency2',8, ...
        %     %     'PassbandRipple1',1, ...
        %     %     'StopbandAttenuation',2, ...
        %     %     'PassbandRipple2',1, ...
        %     %     'SampleRate',meta.samplingRate);
        % 
        % 
        %     lfp = filtfilt(notchFilter50, lfp')';
        %     lfp_mean = mean(lfp,1);
        %     [pxx, freq] = pwelch(lfp_mean, hamming(512), 256, 1024, meta.samplingRate);
        %     figure(f1)
        %     plot(t_lfp(150000:155000), lfp_mean(150000:155000))
        %     figure(f2)
        %     plot(freq, pxx)
        % 
        %     lfp = filtfilt(notchFilter100, lfp')';
        %     lfp_mean = mean(lfp,1);
        %     [pxx, freq] = pwelch(lfp_mean, hamming(512), 256, 1024, meta.samplingRate);
        %     figure(f1);
        %     plot(t_lfp(150000:155000), lfp_mean(150000:155000))
        %     figure(f2);
        %     plot(freq, pxx)
        %     set(gca, 'YScale', 'log')
        % end

        % focus on left and right columns of probe separately, 
        % replace noisy channel data with interpolation, 
        % then smooth across channels and time
        noiseChans = channelLabels > 0;
        for lr = 1:2
            ch = false(size(noiseChans));
            ch(lr:2:end) = true;
            if any(noiseChans(ch))
                [X, T] = ndgrid(find(ch & ~noiseChans), t_lfp);
                F = griddedInterpolant(X, T, ...
                    lfp(ch & ~noiseChans, :), 'linear');
                [X, T] = ndgrid(find(ch & noiseChans), t_lfp);
                lfp(ch & noiseChans,:) = F(X, T);
            end
            lfp(lr:2:end,:) = smoothdata2(lfp(lr:2:end,:), "gaussian", ...
                smoothLFP);
        end
        clear F
        
        %% Determine top of SC base on VEP (visually evoked potential)
        fprintf('  Frame-evoked LFP...\n')
        % Get stimulus triggered LFP
        stimFlat = reshape(stim.frames, size(stim.frames,1), []);
        baselinePerFrame = traces.getAlignedTraces(lfp', t_lfp, ...
            stim.times, baseWin); % [t x frames x channels]
        % [1 x frames x channels]:
        baselinePerFrame = mean(baselinePerFrame, 1, "omitnan"); 
        evokedPerFrame = traces.getAlignedTraces(lfp', t_lfp, ...
            stim.times, stimWin); % [t x frames x channels]
        evokedPerFrame = permute(evokedPerFrame - baselinePerFrame, ...
            [3 1 2]); % [channels x t x frames]
        clear baselinePerFrame
        % average across left and right channel columns
        evokedPerFrame = (evokedPerFrame(1:2:end,:,:) + ...
            evokedPerFrame(2:2:end,:,:)) ./ 2;

        % Get pixel-evoked LFP
        fprintf('  Pixel-evoked LFP...\n')
        lfpEvoked = NaN([size(evokedPerFrame, [1 2]), 2, ...
            size(stimFlat,2)]); % [channels x time x Black/White x pixels]
        for pix = 1:size(stimFlat,2)
            ind_BW = [stimFlat(:,pix) < 0, stimFlat(:,pix) > 0];
            for bw = 1:2
                lfpEvoked(:, :, bw, pix) = ...
                    mean(evokedPerFrame(:,:,ind_BW(:,bw)), 3, "omitnan");
            end
        end

        % Choose response to black or white pixels, whatever is larger
        [~,maxBW] = max(max(abs(lfpEvoked), [], [1 2 4]));
        lfpEvoked = squeeze(lfpEvoked(:,:,maxBW,:)); % [channels x time x pixels]
        % Determine best pixel and time point
        [~, maxPix] = max(max(abs(lfpEvoked), [], [1 2]));
        [~, maxT] = max(max(abs(lfpEvoked), [], [1 3]));
        % Get LFP profile across channels, average across left and right
        % channel columns
        profile = lfpEvoked(:,maxT,maxPix);

        % Find superficial border of SC
        [minAmp, chanMin] = min(profile);
        ind = find(profile(chanMin+1:end) > minAmp * surfaceAmplitude, 1);
        chanTop = ceil(interp1(profile(chanMin:chanMin+ind), chanMin:chanMin+ind, ...
            minAmp * surfaceAmplitude)) * 2;

        %% Determine SGS-SO border based on CSD
        fprintf('  CSD...\n')
        % Calculate CSD
        if maxBW == 1
            indBestPixelFrames = stimFlat(:,maxPix) < 0;
        else
            indBestPixelFrames = stimFlat(:,maxPix) > 0;
        end
        evokedByBestPixel = evokedPerFrame(:,:,indBestPixelFrames);
        clear evokedPerFrame
        csd = mean( diff( diff(evokedByBestPixel,1,1) ,1,1) ,3) ./ (chanDistance^2);
        csd = padarray(csd, 1, 0);

        % locate source and sink, then centre between them
        [~, maxT] = max(max(csd, [], 1));
        [~, chanSource] = max(csd(:,maxT));
        chanSink = chanSource - find(csd(chanSource-1:-1:1, maxT) < 0, 1);
        chan_SGS_SO = floor(interp1(csd(chanSink:chanSource, maxT), ...
            chanSink:chanSource, 0)) * 2;

        %% Save results
        writeNPY([chanTop chan_SGS_SO], fullfile(f, '_ss_recordings.scChannels.npy'))

        %% Make plots
        spikes = io.getEphysData(f);
        
        figure('Position', [3 100 1466 700])
        tiledlayout(1,4)

        nexttile
        m = max(abs(profile),[],"all");
        imagesc(stimWin, channelDepths, lfpEvoked(:,:,maxPix), [-m m])
        colormap(colors)
        c = colorbar;
        c.Label.String ='mV';
        hold on
        plot(stimWin, [1 1].*channelDepths(chanTop/2), 'k', 'LineWidth', 2)
        plot(stimWin, [1 1].*channelDepths(chan_SGS_SO/2), 'k:', 'LineWidth', 2)
        set(gca, 'YDir', 'normal')
        xlabel('Time (ms)')
        ylabel('Distance from probe tip (um)')
        title('Evoked LFP')
        legend('Top of SC', 'SGS-SO-border')

        nexttile
        m = max(abs(csd),[],"all");
        imagesc(stimWin, channelDepths, csd, [-m m])
        colormap(colors)
        c = colorbar;
        c.Label.String = 'mV/mm^2';
        hold on
        plot(stimWin, [1 1].*channelDepths(chanTop/2), 'k', 'LineWidth', 2)
        plot(stimWin, [1 1].*channelDepths(chan_SGS_SO/2), 'k:', 'LineWidth', 2)
        set(gca, 'YDir', 'normal')
        xlabel('Time (ms)')
        title('Current-source density')

        nexttile([1 2])
        ind = spikes.times > stim.times(1) & spikes.times < stim.times(1)+100;
        c = spikes.amps(ind);
        c = c - min(c);
        c = c ./ prctile(c,90);
        c(c>1) = 1;
        scatter(spikes.times(ind) - stim.times(1), ...
            spikes.depths(ind), 1, c, "filled")
        cg = gray;
        colormap(gca, flip(cg))
        hold on
        plot([0 100], [1 1].*channelDepths(chanTop/2), 'r', 'LineWidth', 2)
        plot([0 100], [1 1].*channelDepths(chan_SGS_SO/2), 'r:', 'LineWidth', 2)
        xlim([0 100])
        ylim([min(channelDepths) max(channelDepths)])
        set(gca, 'YDir', 'normal')
        xlabel('Time (ms)')
        title('Spikes')

        sgtitle(sprintf('%s %s: Top: %d, SGS-SO: %d', name, date, ...
            chanTop, chan_SGS_SO))
        saveas(gcf, fullfile(fPlots, sprintf('%s_%s.jpeg', name, date)))
        close gcf
    end
end