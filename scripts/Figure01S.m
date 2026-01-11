function Figure01S(folders, glob)

%% Parameters
% for example traces
buffer = 1; % in sec (before and after stim period)

% for tuning curves
yLims = [-1.1 14.7; -0.8 6.0; -1.6 13.1; -0.4 7.1];

maxP = 0.05;
exp = {'gratingsDrifting', 'bars', 'gratingsStatic'};

%% Examples
% datasets with drifting gratings, static gratings, and bars
ex = {'SS044', '2015-05-15', [176 330 5 364]};

%% For all plots
fPlots = fullfile(folders.plots, 'Figures', 'Figure01S');
if ~isfolder(fPlots)
    mkdir(fPlots)
end

%% Example calcium traces and tuning curves
f = fullfile(folders.data, 'neurons', ex{1,1}, ex{1,2});
% load data
calc = io.getCalciumData(f);
units = ex{3};
for k = 1:3
    stim = io.getGratingInfo(f, exp{k});
    ind = calc.time >= stim.times(1)-1 & calc.time <= stim.times(end)+1;
    % subtract 8th percentile of each trace
    tr = calc.traces(:,units) - prctile(calc.traces(ind,units), 8, 1);

    krnlFits = io.getStimResponseFits(f, exp{k});
    [dirTuning, oriTuning] = io.getTuningResults(f, exp{k});
    % stim. duration and non-gray stimuli
    stimDur = median(diff(stim.times,1,2));
    if k < 3
        stimPars = stim.directions;
    else
        stimPars = stim.orientations;
    end
    validStims = find(~isnan(stimPars));

    % align calcium trace to stimuli
    [alignedTrace, t] = traces.getAlignedTraces( ...
        tr, calc.time, ...
        stim.times(:,1), [-buffer stimDur+buffer]);
    % align predicted trace to stimuli
    alignedPredicted = traces.getAlignedTraces( ...
        krnlFits.prediction(:,units), krnlFits.time_prediction, ...
        stim.times(:,1), [-buffer stimDur+buffer]);
    % loop over all examples
    for iUnit = 1:length(units)
        mini = min(alignedTrace(:,:,iUnit), [], "all");
        maxi = max(alignedTrace(:,:,iUnit), [], "all");
        figure('Position',[3 570 1915 420]);
        tiledlayout(1, length(validStims), "TileSpacing", "tight", ...
            "Padding", "tight")
        % loop over all stimuli
        for st = 1:length(validStims)
            indSt = stim.ids == validStims(st);
            nexttile
            hold on
            % mark time of stimulus
            fill([0 0 stimDur stimDur], [mini maxi maxi mini], 'k', ...
                'FaceColor', [1 1 1].*0.9, 'EdgeColor', 'none')
            % single trial traces
            plot(t, alignedTrace(:,indSt,iUnit), ...
                "Color", [1 1 1].*0.5);
            % mean predicted traces (from kernel fit)
            plot(t, mean(alignedPredicted(:,indSt,iUnit), 2, "omitnan"), ...
                'k', "LineWidth", 1)
            set(gca, "Box", "off")
            xlim(t([1 end]))
            ylim([mini maxi])
            if st == 1
                xlabel('Time (s)')
                ylabel('\DeltaF/F')
            end
            if k < 3
                title(sprintf('%d deg', stimPars(validStims(st))))
            else
                title(sprintf('%d deg (phase: %d deg)', ...
                    stimPars(validStims(st)), stim.phases(validStims(st))))
            end
        end
        sgtitle(sprintf('ROI %03d', units(iUnit)))
        io.saveFigure(gcf, fPlots, sprintf('example_%s_stimTraces_%s_%s_%03d', ...
            exp{k}, ex{1}, ex{2}, units(iUnit)))

        % plot direction tuning curve
        if k < 3
            amps = krnlFits.amplitudes(:, validStims([1:end 1]), units(iUnit));
            drct = [stimPars(validStims); 360]';
            figure('Position', [700 150 475 325]);
            hold on
            plot(drct + randn(size(amps)).*3, ...
                amps, '.', 'Color', [1 1 1].*0.5, 'MarkerSize', 10);
            plot(drct, mean(amps,1,'omitnan'), ...
                'k.-', 'MarkerSize', 30, 'LineWidth', 1);
            if dirTuning.pValue(units(iUnit)) < maxP
                plot(dirTuning.preference(units(iUnit)), maxi, 'v', 'MarkerSize', 8, ...
                    'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'none');
            end
            set(gca, 'box', 'off', 'XTick', 0:90:360)
            xlim([-10 370])
            ylim(yLims(k,:))
            xlabel('Direction (deg)')
            ylabel('\DeltaF/F (kernel amplitude)')
            title(sprintf('ROI %03d: Direction tuning', units(iUnit)))
            io.saveFigure(gcf, fPlots, ...
                sprintf('example_%s_dirTuningCurve_%s_%s_%03d', ...
                exp{k}, ex{1}, ex{2}, units(iUnit)))
        end

        % plot orientation tuning curve
        amps = krnlFits.amplitudes(:, validStims, units(iUnit));
        if k < 3
            amps = reshape(amps, size(amps,1), size(amps,2)/2, 2);
            amps = reshape(permute(amps, [1 3 2]), size(amps,1)*2, []);
            ortn = [stimPars(validStims(1:length(validStims)/2)); 180]';
        else
            stimUni = unique(stimPars(validStims));
            tmp = [];
            for j = 1:length(stimUni)
                tmp = [tmp, reshape(amps(:,stimPars==stimUni(j)),[],1)];
            end
            amps = tmp;
            ortn = [stimUni; 180]';
        end
        amps = amps(:,[1:end 1]);
        figure('Position', [1250 150 475 325]);
        hold on
        plot(ortn + randn(size(amps)).*3, ...
            amps, '.', 'Color', [1 1 1].*0.5, 'MarkerSize', 10);
        plot(ortn, mean(amps,1,'omitnan'), ...
            'k.-', 'MarkerSize', 30, 'LineWidth', 1);
        if oriTuning.pValue(units(iUnit)) < maxP
            plot(oriTuning.preference(units(iUnit)), maxi, 'v', 'MarkerSize', 8, ...
                'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'none');
        end
        set(gca, 'box', 'off', 'XTick', 0:45:180)
        xlim([-10 190])
        ylim(yLims(k,:))
        xlabel('Orientation (deg)')
        ylabel('\DeltaF/F (kernel amplitude)')
        title(sprintf('ROI %03d: Orientation tuning', units(iUnit)))
        io.saveFigure(gcf, fPlots, sprintf('example_%s_oriTuningCurve_%s_%s_%03d', ...
            exp{k}, ex{1}, ex{2}, units(iUnit)))
    end
end

%% Preferences across stimulus paradigms
% Collect data
subjDirs = dir(fullfile(folders.data, 'neurons', 'SS*'));
R2 = [];
P = [];
dirPreferences = [];
oriPreferences = [];
dirTuned = [];
oriTuned = [];
nRecorded = [0 0 0];
nResponsive = [0 0 0];
nDirTuned = [0 0 0];
nOriTuned = [0 0 0];
nSessions = [0 0 0];
for subj = 1:length(subjDirs) % animals
    name = subjDirs(subj).name;
    dateDirs = dir(fullfile(folders.data, 'neurons', name, '2*'));
    for d = 1:length(dateDirs) %dates
        date = dateDirs(d).name;
        f = fullfile(folders.data, 'neurons', name, date);
        r2 = [];
        pKernel = [];
        for k = 1:3
            % ignore session if stimulus was not presented
            if ~isfile(fullfile(f, sprintf('_ss_%s.intervals.npy', exp{k})))
                continue
            end

            krnlFits = io.getStimResponseFits(f, exp{k});
            [dirTuning, oriTuning] = io.getTuningResults(f, exp{k});
            numUnits = length(krnlFits.pValue);
            if isempty(r2)
                r2 = NaN(numUnits,3);
                pKernel = NaN(numUnits,3);
                dp = NaN(numUnits,2);
                op = NaN(numUnits,3);
                dt = NaN(numUnits,2);
                ot = NaN(numUnits,3);
            end

            r2(:,k) = krnlFits.R2;
            pKernel(:,k) = krnlFits.pValue;
            if ~isempty(dirTuning)
                dp(:,k) = dirTuning.preference;
                dt(:,k) = dirTuning.pValue < maxP;
                nDirTuned(k) = nDirTuned(k) + sum(dirTuning.pValue < maxP);
            end
            op(:,k) = oriTuning.preference;
            ot(:,k) = oriTuning.pValue < maxP;

            nRecorded(k) = nRecorded(k) + numUnits;
            nResponsive(k) = nResponsive(k) + ...
                sum(krnlFits.pValue < maxP & krnlFits.R2 > 0);
            nOriTuned(k) = nOriTuned(k) + sum(oriTuning.pValue < maxP);
            nSessions(k) = nSessions(k) + 1;
        end
        R2 = [R2; r2];
        P = [P; pKernel];
        dirPreferences = [dirPreferences; dp];
        oriPreferences = [oriPreferences; op];
        dirTuned = [dirTuned; dt];
        oriTuned = [oriTuned; ot];

        if strcmp(name, ex{1}) && strcmp(date, ex{2})
            n = length(R2) - numUnits;
            indExamples = n + ex{3};
        end
    end
end

for k = 1:3
    fprintf('%s:\n', exp{k})
    fprintf('  %d neurons were recorded across %d sessions\n', ...
        nRecorded(k), nSessions(k))
    ind = P(:,k) < maxP & R2(:,k) > 0;
    fprintf('  Explained variance of kernel-based predictions (mean +- STD): %.4f +- %.4f\n', ...
        mean(R2(ind, k)), std(R2(ind, k)))
    fprintf('  Of %d responsive neurons, %d were tuned to direction, %d were tuned to orientation\n', ...
        nResponsive(k), nDirTuned(k), nOriTuned(k))
end

dotCols = lines(length(ex{3}));

% Scatterplot: preferred direction from drifting gratings vs bars
ind = all(dirTuned(:, [1 2]) == 1, 2);
diffs = dirPreferences(ind,1) - dirPreferences(ind,2);
diffs(diffs < -180) = 360 + diffs(diffs < -180);
diffs(diffs > 180) = diffs(diffs > 180) - 360;
m = mean(diffs);
p = signrank(diffs);
% [~,p] = ttest(diffs);
figure('Position', glob.figPositionDefault)
hold on
plot([0 360], [0 360], 'Color', [1 1 1].*0.5)
scatter(dirPreferences(ind,1), dirPreferences(ind,2), 15, 'k', 'filled')
scatter(dirPreferences(indExamples,1), dirPreferences(indExamples,2), ...
    40, dotCols, 'filled')
axis equal
xlim([-10 370])
ylim([-10 370])
set(gca, "Box", "off", "XTick", 0:90:360, "YTick", 0:90:360)
xlabel(exp{1})
ylabel(exp{2})
title(sprintf('Preferred directions (n=%d, <x-y>=%.1f, p=%.3f)', ...
    sum(ind), m, p))
io.saveFigure(gcf, fPlots, sprintf('prefDir_%s-%s', exp{1}, exp{2}));

% Scatterplot: preferred orientation from drifting gratings vs bars
ind = all(oriTuned(:, [1 2]) == 1, 2);
diffs = oriPreferences(ind,1) - oriPreferences(ind,2);
diffs(diffs < -90) = 180 + diffs(diffs < -90);
diffs(diffs > 90) = diffs(diffs > 90) - 180;
m = mean(diffs);
p = signrank(diffs);
% [~,p] = ttest(diffs);
figure('Position', glob.figPositionDefault)
hold on
plot([0 360], [0 360], 'Color', [1 1 1].*0.5)
scatter(oriPreferences(ind,1), oriPreferences(ind,2), 15, 'k', 'filled')
scatter(oriPreferences(indExamples,1), oriPreferences(indExamples,2), ...
    40, dotCols, 'filled')
axis equal
xlim([-10 190])
ylim([-10 190])
set(gca, "Box", "off", "XTick", 0:45:180, "YTick", 0:45:180)
xlabel(exp{1})
ylabel(exp{2})
title(sprintf('Preferred orientations (n=%d, <x-y>=%.1f, p=%.3f)', ...
    sum(ind), m, p))
io.saveFigure(gcf, fPlots, sprintf('prefOri_%s-%s', exp{1}, exp{2}));

% Scatterplot: preferred orientation from drifting vs static gratings
ind = all(oriTuned(:, [1 3]) == 1, 2);
diffs = oriPreferences(ind,1) - oriPreferences(ind,3);
diffs(diffs < -90) = 180 + diffs(diffs < -90);
diffs(diffs > 90) = diffs(diffs > 90) - 180;
m = mean(diffs);
p = signrank(diffs);
% [~,p] = ttest(diffs);
figure('Position', glob.figPositionDefault)
hold on
plot([0 360], [0 360], 'Color', [1 1 1].*0.5)
scatter(oriPreferences(ind,1), oriPreferences(ind,3), 15, 'k', 'filled')
scatter(oriPreferences(indExamples,1), oriPreferences(indExamples,3), ...
    40, dotCols, 'filled')
axis equal
xlim([-10 190])
ylim([-10 190])
set(gca, "Box", "off", "XTick", 0:45:180, "YTick", 0:45:180)
xlabel(exp{1})
ylabel(exp{3})
title(sprintf('Preferred orientations (n=%d), <x-y>=%.1f, p=%.3f', ...
    sum(ind), m, p))
io.saveFigure(gcf, fPlots, sprintf('prefOri_%s-%s', exp{1}, exp{3}));

% Scatterplot: preferred orientation from static gratings vs bars
ind = all(oriTuned(:, [3 2]) == 1, 2);
diffs = oriPreferences(ind,3) - oriPreferences(ind,2);
diffs(diffs < -90) = 180 + diffs(diffs < -90);
diffs(diffs > 90) = diffs(diffs > 90) - 180;
m = mean(diffs);
p = signrank(diffs);
% [~,p] = ttest(diffs);
figure('Position', glob.figPositionDefault)
hold on
plot([0 360], [0 360], 'Color', [1 1 1].*0.5)
scatter(oriPreferences(ind,3), oriPreferences(ind,2), 15, 'k', 'filled')
scatter(oriPreferences(indExamples,3), oriPreferences(indExamples,2), ...
    40, dotCols, 'filled')
axis equal
xlim([-10 190])
ylim([-10 190])
set(gca, "Box", "off", "XTick", 0:45:180, "YTick", 0:45:180)
xlabel(exp{3})
ylabel(exp{2})
title(sprintf('Preferred orientations (n=%d), <x-y>=%.1f, p=%.3f', ...
    sum(ind), m, p))
io.saveFigure(gcf, fPlots, sprintf('prefOri_%s-%s', exp{3}, exp{2}));