function Figure06_examplesTuning(folders, fPlots, animal, date, units)

buffer = 1; % in sec (before and after stim period)
baselineTime = 0.5;
stimulusDur = 2;
durationSlack = 0.01;
binSize = 0.02;
sigma = 0.1;

f = fullfile(folders.data, 'ephys', animal, date);

spikeData = io.getEphysData(f);
[dirTuning, oriTuning] = io.getTuningResults(f, 'gratingsDrifting');
stimData = io.getGratingInfo(f, 'gratingsDrifting');
stimData = stimuli.fixNonDirectionFeatures(stimData);
validStims = find(~isnan(stimData.directions));
stimReps = max(histcounts(stimData.directions(stimData.ids), ...
    [stimData.directions(validStims); 361]));
stimDur = diff(stimData.times, 1, 2);
% identify trials that are too short
invalidTrials = find(stimDur < stimulusDur - durationSlack);
stimDur = median(stimDur);

for k = 1:length(units)
    indUnit = find(spikeData.clusterIDs == units(k));
    t = spikeData.times(spikeData.clusters == units(k));
    [sp_aligned, trial] = events.alignData(t, ...
        stimData.times(:,1), [-buffer, stimDur + buffer]);
    invalid = ismember(trial, invalidTrials);
    sp_aligned(invalid) = [];
    trial(invalid) = [];
    [traces_singleTrials, bins] = events.tracesFromEvents(sp_aligned, trial, ...
        (1:length(stimData.ids))', [-buffer, stimDur + buffer], ...
        binSize, sigma, false);

    mini = min(traces_singleTrials, [], "all");
    maxi = max(traces_singleTrials, [], "all");
    figure('Position',[3 570 1915 420]);
    tiledlayout(1, length(validStims), "TileSpacing", "tight", ...
        "Padding", "tight")
    amplitudes = NaN(stimReps, length(validStims));
    % loop over all stimuli
    for st = 1:length(validStims)
        stimTrials = find(stimData.ids == validStims(st));
        nexttile
        hold on
        % mark time of stimulus
        fill([0 0 stimDur stimDur], [mini maxi maxi mini], 'k', ...
            'FaceColor', [1 1 1].*0.9, 'EdgeColor', 'none')
        % single trial traces
        plot(bins, traces_singleTrials(stimTrials,:), "Color", [1 1 1].*0.5);
        % mean trace
        plot(bins, mean(traces_singleTrials(stimTrials,:), 1, "omitnan"), ...
            'k', "LineWidth", 1)
        set(gca, "Box", "off")
        xlim(bins([1 end]))
        ylim([mini maxi])
        if st == 1
            xlabel('Time (s)')
            ylabel('Firing rate (sp/s)')
        end
        title(sprintf('%d deg', stimData.directions(validStims(st))))
        
        for rep = 1:length(stimTrials)
            tr = stimTrials(rep);
            if ismember(tr, invalidTrials)
                continue
            end
            ind_spikes = trial == tr;
            amplitudes(rep, st) = ...
                sum(sp_aligned(ind_spikes) >= 0 & ...
                sp_aligned(ind_spikes) < stimDur) / ...
                stimDur - ...
                sum(sp_aligned(ind_spikes) >= -baselineTime & ...
                sp_aligned(ind_spikes) < 0) / ...
                baselineTime;
        end
    end
    sgtitle(sprintf('Unit %04d: DS = %.2f (p=%.3f), OS = %.2f (p=%.3f)', ...
        units(k), dirTuning.selectivity(indUnit), ...
        dirTuning.pValue(indUnit), oriTuning.selectivity(indUnit), ...
        oriTuning.pValue(indUnit)))
    io.saveFigure(gcf, fPlots, sprintf('example_stimTraces_%s_%s_%04d', ...
        animal, date, units(k)))

    % plot direction tuning curve
    amps = amplitudes(:, [1:end 1]);
    drct = [stimData.directions(validStims); 360]';
    mini = min(amps, [], "all");
    maxi = max(amps, [], "all");
    range = maxi - mini;
    figure('Position', [700 150 475 325])
    hold on
    plot(drct + randn(size(amps)).*3, ...
        amps, '.', 'Color', [1 1 1].*0.5, 'MarkerSize', 10);
    plot(drct, mean(amps,1,'omitnan'), ...
        'k.-', 'MarkerSize', 30, 'LineWidth', 1);
    plot(dirTuning.preference(indUnit), maxi, 'v', 'MarkerSize', 8, ...
        'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'none');
    set(gca, 'box', 'off', 'XTick', drct(1:3:end))
    xlim([-10 370])
    ylim([min([0 mini-0.05*range]) maxi+0.05*range])
    xlabel('Direction (deg)')
    ylabel('Firing rate (sp/s)')
    title(sprintf('Unit %04d: Direction tuning', units(k)))
    io.saveFigure(gcf, fPlots, sprintf('example_dirTuningCurve_%s_%s_%04d', ...
        animal, date, units(k)))

    % plot orientation tuning curve
    amps = amplitudes;
    amps = reshape(amps, size(amps,1), size(amps,2)/2, 2);
    amps = reshape(permute(amps, [1 3 2]), size(amps,1)*2, []);
    amps = amps(:,[1:end 1]);
    ortn = [stimData.directions(validStims(1:length(validStims)/2)); 180]';
    mini = min(amps, [], "all");
    maxi = max(amps, [], "all");
    range = maxi - mini;
    figure('Position', [1250 150 475 325])
    hold on
    plot(ortn + randn(size(amps)).*3, ...
        amps, '.', 'Color', [1 1 1].*0.5, 'MarkerSize', 10);
    plot(ortn, mean(amps,1,'omitnan'), ...
        'k.-', 'MarkerSize', 30, 'LineWidth', 1);
    plot(oriTuning.preference(indUnit), maxi, 'v', 'MarkerSize', 8, ...
        'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'none');
    set(gca, 'box', 'off', 'XTick', ortn(1:3:end))
    xlim([-10 190])
    ylim([min([0 mini-0.05*range]) maxi+0.05*range])
    xlabel('Orientation (deg)')
    ylabel('Firing rate (sp/s)')
    title(sprintf('Unit %04d: Orientation tuning', units(k)))
    io.saveFigure(gcf, fPlots, sprintf('example_oriTuningCurve_%s_%s_%04d', ...
        animal, date, units(k)))
end