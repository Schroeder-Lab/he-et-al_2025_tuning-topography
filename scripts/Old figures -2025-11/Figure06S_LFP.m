function Figure06S_LFP(folders, fPlots, animal, date)

%% Parameters
smoothLFP = {5 21};
stimWin = [-.05 .15];
baseWin = [-.05 0];
surfaceAmplitude = 0.25;
colors = colmaps.getBlueWhiteRedMap(255);
window = [100 200];

%% Load data
f = fullfile(folders.data, 'ephys', animal, date);
% Note: raw data are not publicly shared due to large file size
f_raw = fullfile(folders.dataRawEphys, animal, date, 'exp01', ...
    sprintf('%s_%s_g0', animal, date), ...
    sprintf('%s_%s_g0_imec0', animal, date));

% load visual stimulus data
stim = io.getVisNoiseInfo(f);

% load LFP meta data
meta = io.getLFPMetaData(fullfile(f_raw, ...
    sprintf('%s_%s_g0_t0.imec0.lf.meta', animal, date)));
channelLabels = readNPY(fullfile(f_raw, 'ibl_sorter', ...
    'channel_labels.npy'));
channelCoords = readNPY(fullfile(f, 'channels.localCoordinates.npy'));
channelDepths = channelCoords(1:2:end,2);
chanDistance = diff(channelDepths(1:2)) / 1000; % millimeters

% load LFP data
fileLFP = dir(fullfile(f_raw, ...
    sprintf('%s_%s_g0_t0.imec0.lf.bin', animal, date)));
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

% subtract mean lfp trace (across channels) from each channel
lfp_mean = mean(lfp, 1);
lfp = lfp - lfp_mean;

%% Determine visually evoked potential
% used to determine top of SC
% Get stimulus-aligned LFP, subtract baseline
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
stimFlat = reshape(stim.frames, size(stim.frames,1), []);
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

%% Determine CSD
% used to determine SGS-SO border
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
xlabel('Time (s)')
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
xlabel('Time (s)')
title('Current-source density')

nexttile([1 2])
ind = spikes.times > stim.times(1)+window(1) & ...
    spikes.times < stim.times(1)+window(2);
c = spikes.amps(ind);
c = c - min(c);
c = c ./ prctile(c,90);
c(c>1) = 1;
scatter(spikes.times(ind) - stim.times(1), ...
    spikes.depths(ind), 1, c, "filled")
cg = gray;
colormap(gca, flip(cg))
hold on
plot(window, [1 1].*channelDepths(chanTop/2), 'r', 'LineWidth', 2)
plot(window, [1 1].*channelDepths(chan_SGS_SO/2), 'r:', 'LineWidth', 2)
xlim(window)
ylim([min(channelDepths) max(channelDepths)])
set(gca, 'YDir', 'normal')
xlabel('Time (s)')
title('Spikes')

sgtitle(sprintf('%s %s: Top: %d, SGS-SO: %d', animal, date, ...
    chanTop, chan_SGS_SO))
io.saveFigure(gcf, fPlots, sprintf('LFP_%s_%s', animal, date))