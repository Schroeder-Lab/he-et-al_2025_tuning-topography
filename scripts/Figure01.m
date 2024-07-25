%% Folders
getFolders;

%% Parameters

%% Examples
ex = cell(2,4); % rows: (1) bouton, (2) neuron, columns: (1) ori, (2) dir
ex(1,1,:) = {'SS076', '2017-10-02', 6, 4};
ex(1,2,:) = {'SS076', '2017-10-02', 6, 6};
ex(2,:) = {'SS044', '2015-04-28', 3, [227 231 236 240 244 245 249 252 256 258 260 262 263 264 268 273 281 285 293 313 321 328 369 378 379 380 383 385 393 393 402 408 412 425 426 430 440]};
sets = {'boutons', 'neurons'};

%% Add paths
addpath(genpath(fullfile(folders.tools, 'npy-matlab')))
addpath(fullfile(folders.repo))

%% Example calcium traces and tuning curves
buffer = 1; % in sec (before and after stim period)
for s = 1:2
    str = sets{k};
    f = fullfile(folders.data, str, ex{s,1}, ex{s,2});
    calc = io.getCalciumData(f);
    stim = io.getGratingInfo(f, 'gratingsDrifting');
    krnlFits = io.getStimResponseFits(f, 'gratingsDrifting');

    % loop over all examples
    units = ex{s,4};
    for k = 1:length(units)
        % align calcium traces to stimuli

        % align prediction traces to stimuli

        % plot for each stimulus:
        % single-trial traces +
        % mean of prediction
    end

    cellID = examplesTun{ex,3}==planes & examplesTun{ex,4}==ids;
    t = time + planeDelays(examplesTun{ex,3});
    stimMatrix = exp.buildStimMatrix(stimSequence, stimIntervals, t);
    directions(isnan(directions)) = [];
    timeBin = median(diff(time));
    repetitions = sum(stimSequence == 1);
    stimDurInFrames = round(sum(stimMatrix(1,:)) / repetitions);
    stimDur = stimDurInFrames * timeBin;
    offset = ceil(buffer / timeBin);
    resp = squeeze(exp.getTracesPerStimulus(traces(:,cellID), stimMatrix, ...
        [1 1] .* offset)); % [stimulus x trial x time]
    predResp = squeeze(exp.getTracesPerStimulus(predictions(:,cellID), ...
        stimMatrix, [1 1] .* offset));

    ind = ~largePupil;
    ind = repmat(ind, 1, 1, size(resp,3));
    temp = resp(1:end-1,:,:);
    temp(~ind) = NaN;
    respMean = nanmean(temp, 2);
    temp = predResp(1:end-1,:,:);
    temp(~ind) = NaN;
    predMean = nanmean(temp, 2);
    respMean = respMean(1:3:end,:);
    predMean = predMean(1:3:end,:);

    mini = min([respMean(:); predMean(:)]);
    maxi = max([respMean(:); predMean(:)]);
    rng = maxi - mini;
    mini = mini - 0.05*rng;
    maxi = maxi + 0.05*rng;
    xDist = .5;
    traceDur = stimDur + 2*buffer;
    respTime = (-offset:stimDurInFrames+offset-1) .* timeBin;

    figure('Position',[50 500 1000 420])
    hold on
    h = [0 0];
    x0 = 0;
    for st = 1:size(respMean,1)
        fill([0 stimDur stimDur 0] + x0, ...
            [mini mini maxi maxi], 'k', 'FaceColor', 'k', ...
            'FaceAlpha', 0.1, 'EdgeColor', 'none')
        plot(respTime([1 end]) + x0, [0 0], 'k:')
        h1 = plot(respTime + x0, squeeze(respMean(st,:)), ...
            'Color', 'k');
        h2 = plot(respTime + x0, squeeze(predMean(st,:)), ...
            'Color', 'k', 'LineWidth', 2);
        x0 = x0 + traceDur + xDist;
        if st==1
            h(1) = h2;
            h(2) = h1;
        end
    end
    axis tight
    set(gca, 'XTick', [0 stimDur])
    legend(h, {'Prediction from kernel fit','Data'})
    xlabel('Stimuli')
    ylabel('\DeltaF/F')
    title(sprintf('Example %d',ex))

    figure
    hold on
    amps = amplitudes(:,:,cellID);
    amps(largePupil) = NaN;
    plot(1:360, curves(cellID,:), 'Color', 'k', 'LineWidth',2);
    m = nanmean(amps,2);
    s = nanstd(amps,0,2) ./ sqrt(sum(~largePupil,2));
    errorbar([directions; 360], m([1:end 1]), s([1:end 1]), 'o', ...
        'Color', 'k', 'CapSize', 2, 'MarkerFaceColor', 'k')
    plot([0 360], [0 0], 'k:', 'LineWidth', 2)
    set(gca,'XTick',0:90:360)
    title(sprintf('Example %d',ex))
    xlim([0 360])
    mini = min([mini; m-s; curves(cellID,:)']);
    maxi = max([maxi; m+s; curves(cellID,:)']);
    ylim([mini maxi])
    xlabel('Direction (in degrees)')
    ylabel('\DeltaF/F')
end

%% Population direction tuning curves

%% Population orientation tuning curves

%% Total direction preferences histogram

%% Total orientation preferences histogram

%% Direction versus orientation selectivity indices