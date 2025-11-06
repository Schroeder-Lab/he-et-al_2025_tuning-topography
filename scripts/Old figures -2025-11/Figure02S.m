function Figure02S(folders)

%% Parameters
sets = {'boutons', 'neurons'};
selectivityThresholds = [0.2 0.2; 0.1 0.2];
minROIs = 15;
binSize = [5, 20];
stepSize = [2.5, 5];
xLims = [50 500];

exp = {'bars', 'gratingsStatic'};
maxP = 0.05;
numPerm = 1000;

%% For all plots
fPlots = fullfile(folders.plots, 'Figures', 'Figure02S');
if ~isfolder(fPlots)
    mkdir(fPlots)
end

%% Plot pairwise distance versus tuning difference : DS- and OS-only units
for s = 1:2 % boutons and neurons
    subjDirs = dir(fullfile(folders.data, sets{s}, 'SS*'));
    dirDist = {};
    dirDiff = {};
    dirDiffNull = {};
    oriDist = {};
    oriDiff = {};
    oriDiffNull = {};
    rec = 1;
    for subj = 1:length(subjDirs) % animals
        name = subjDirs(subj).name;
        fprintf('%s\n', name)
        dateDirs = dir(fullfile(folders.data, sets{s}, name, '2*'));
        for dt = 1:length(dateDirs) %dates
            date = dateDirs(dt).name;
            f = fullfile(folders.data, sets{s}, name, date);
            % ignore session if stimulus was not presented
            if ~isfile(fullfile(f, '_ss_gratingsDrifting.intervals.npy'))
                continue
            end
                
            % load data
            data = io.getRecordingInfo(f);
            roiPos = data.roiPositions(:,1:2);
            [dirTuning, oriTuning] = io.getTuningResults(f, 'gratingsDrifting');

            dp = dirTuning.preference;
            validDir = ~isnan(dp) & dirTuning.pValue < maxP & ...
                dirTuning.selectivity >= selectivityThresholds(1,1) & ...
                oriTuning.selectivity <= selectivityThresholds(1,2);
            if sum(validDir) < minROIs
                dirDist{rec} = [];
                dirDiff{rec} = [];
                dirDiffNull{rec} = [];
            else
                % for all unit pairs, determine distance in brain (ignore
                % depth);
                ddist = spatial.determineDistance(roiPos(validDir,1), ...
                    roiPos(validDir,2));
                dp = dp(validDir);
                % for all unit pairs, determine difference between preferred
                % directions
                ddiff = tuning.determinePreferenceDiff(dp, 'dir');
                % permute preferences to test significance
                ddiffPermuted = NaN(length(ddiff), numPerm);
                rng('default');
                for k = 1:numPerm
                    order = randperm(length(dp));
                    ddiffPermuted(:,k) = tuning.determinePreferenceDiff(dp(order), 'dir');
                end
                % collect results
                dirDist{rec} = ddist;
                dirDiff{rec} = ddiff;
                dirDiffNull{rec} = ddiffPermuted;
            end
            op = oriTuning.preference;
            validOri = ~isnan(op) & oriTuning.pValue < maxP & ...
                oriTuning.selectivity >= selectivityThresholds(2,1) & ...
                dirTuning.selectivity <= selectivityThresholds(2,2);
            if sum(validOri) < minROIs
                oriDist{rec} = [];
                oriDiff{rec} = [];
                oriDiffNull{rec} = [];
            else
                % for all unit pairs, determine distance in brain (ignore
                % depth);
                odist = spatial.determineDistance(roiPos(validOri,1), ...
                    roiPos(validOri,2));
                op = op(validOri);
                % for all unit pairs, determine difference between preferred
                % orientations
                odiff = tuning.determinePreferenceDiff(op, 'ori');
                % permute preferences to test significance
                odiffPermuted = NaN(length(odiff), numPerm);
                rng('default');
                for k = 1:numPerm
                    order = randperm(length(op));
                    odiffPermuted(:,k) = tuning.determinePreferenceDiff(op(order), 'ori');
                end
                % collect results
                oriDist{rec} = odist;
                oriDiff{rec} = odiff;
                oriDiffNull{rec} = odiffPermuted;
            end
            rec = rec + 1;
        end
    end

    % plot distance vs tuning difference across all datasets
    n = sum(~any(isnan([cat(1, dirDist{:}) cat(1, dirDiff{:})]), 2));
    fig = spatial.plotPrefDiffVsDist(cat(1, dirDist{:}), ...
        cat(1, dirDiff{:}), cat(1, dirDiffNull{:}), ...
        binSize(s), stepSize(s), false);
    set(gca, 'YTick', 0:45:180)
    xlim([0 xLims(s)])
    ylim([0 180])
    title(['\DeltaDirection pref. vs \Deltaposition (n = ' num2str(n) ')'])
    io.saveFigure(fig, fPlots, ...
        sprintf('distanceAll_%s_directionOnly', sets{s}))
    n = sum(~any(isnan([cat(1, oriDist{:}) cat(1, oriDiff{:})]), 2));
    fig = spatial.plotPrefDiffVsDist(cat(1, oriDist{:}), ...
        cat(1, oriDiff{:}), cat(1, oriDiffNull{:}), ...
        binSize(s), stepSize(s), false);
    set(gca, 'YTick', 0:45:90)
    xlim([0 xLims(s)])
    ylim([0 90])
    title(['\DeltaOrientation pref. vs \Deltaposition (n = ' num2str(n) ')'])
    io.saveFigure(fig, fPlots, ...
        sprintf('distanceAll_%s_orientationOnly', sets{s}))
end

%% Plot pairwise distance in brain versus difference in tuning preference 
% measured in response to bars and static gratings

subjDirs = dir(fullfile(folders.data, 'neurons', 'SS*'));
dist = {};
dirDiff_bars = {};
dirDiffNull_bars = {};
oriDiff_bars = {};
oriDiffNull_bars = {};
oriDiff_static = {};
oriDiffNull_static = {};
rec = 1;
for subj = 1:length(subjDirs) % animals
    name = subjDirs(subj).name;
    dateDirs = dir(fullfile(folders.data, 'neurons', name, '2*'));
    for dt = 1:length(dateDirs) % dates
        date = dateDirs(dt).name;
        f = fullfile(folders.data, 'neurons', name, date);

        % ignore session if neither stimulus was presented
        if ~isfile(fullfile(f, sprintf('_ss_%s.intervals.npy', exp{1}))) && ...
                ~isfile(fullfile(f, sprintf('_ss_%s.intervals.npy', exp{2})))
            continue
        end

        % load data
        data = io.getRecordingInfo(f);
        roiPos = data.roiPositions(:,1:2);
        % for all unit pairs, determine distance in brain (ignore depth);
        dist{rec} = spatial.determineDistance(roiPos(:,1), roiPos(:,2));
        for k = 1:2
            % ignore session if stimulus was not presented
            if ~isfile(fullfile(f, sprintf('_ss_%s.intervals.npy', exp{k})))
                if k == 1 % bars
                    dirDiff_bars{rec} = [];
                    dirDiffNull_bars{rec} = [];
                    oriDiff_bars{rec} = [];
                    oriDiffNull_bars{rec} = [];
                else % static
                    oriDiff_static{rec} = [];
                    oriDiffNull_static{rec} = [];
                end
                continue
            end
            [dirTuning, oriTuning] = io.getTuningResults(f, exp{k});
            if k == 1 % bars
                dp = dirTuning.preference;
                invalid = dirTuning.pValue >= maxP;
                dp(invalid) = NaN;
                % for all unit pairs, determine difference between preferred
                % directions
                ddiff = tuning.determinePreferenceDiff(dp, 'dir');
                % permute preferences to test significance
                ddiffPermuted = NaN(length(ddiff), numPerm);
                rng('default');
                for n = 1:numPerm
                    order = randperm(length(dp));
                    ddiffPermuted(:,n) = ...
                        tuning.determinePreferenceDiff(dp(order), 'dir');
                end
                % collect results
                dirDiff_bars{rec} = ddiff;
                dirDiffNull_bars{rec} = ddiffPermuted;
            end

            op = oriTuning.preference;
            invalid = oriTuning.pValue >= maxP;
            op(invalid) = NaN;
            % for all unit pairs, determine difference between preferred
            % orientations
            odiff = tuning.determinePreferenceDiff(op, 'ori');
            % permute preferences to test significance
            odiffPermuted = NaN(length(odiff), numPerm);
            rng('default');
            for n = 1:numPerm
                order = randperm(length(op));
                odiffPermuted(:,n) = ...
                    tuning.determinePreferenceDiff(op(order), 'ori');
            end
            % collect results
            if k == 1 % bars
                oriDiff_bars{rec} = odiff;
                oriDiffNull_bars{rec} = odiffPermuted;
            else % static
                oriDiff_static{rec} = odiff;
                oriDiffNull_static{rec} = odiffPermuted;
            end
        end
        rec = rec + 1;
    end
end

% Tuning to bars
% 1. plot distance vs direction difference
ind = ~cellfun(@isempty, dirDiff_bars);
n = sum(~any(isnan([cat(1, dist{ind}) cat(1, dirDiff_bars{ind})]), 2));
spatial.plotPrefDiffVsDist(cat(1, dist{ind}), ...
    cat(1, dirDiff_bars{ind}), cat(1, dirDiffNull_bars{ind}), ...
    binSize(2), stepSize(2), false);
xlim([0 400])
ylim([0 180])
clim([0 3.5e-5])
title(['\DeltaDirection pref.: ' exp{1} ' (n = ' num2str(n) ')'])
io.saveFigure(gcf, fPlots, sprintf('dirPref_%s', exp{1}))

% 2. plot distance vs orientation difference
ind = ~cellfun(@isempty, oriDiff_bars);
n = sum(~any(isnan([cat(1, dist{ind}) cat(1, oriDiff_bars{ind})]), 2));
spatial.plotPrefDiffVsDist(cat(1, dist{ind}), ...
    cat(1, oriDiff_bars{ind}), cat(1, oriDiffNull_bars{ind}), ...
    binSize(2), stepSize(2), false);
xlim([0 400])
ylim([0 90])
clim([0 7e-5])
title(['\DeltaOrientation pref.: ' exp{1} ' (n = ' num2str(n) ')'])
io.saveFigure(gcf, fPlots, sprintf('oriPref_%s', exp{1}))

% Tuning to static gratings
% plot distance vs orientation difference
ind = ~cellfun(@isempty, oriDiff_static);
n = sum(~any(isnan([cat(1, dist{ind}) cat(1, oriDiff_static{ind})]), 2));
spatial.plotPrefDiffVsDist(cat(1, dist{ind}), ...
    cat(1, oriDiff_static{ind}), cat(1, oriDiffNull_static{ind}), ...
    binSize(2), stepSize(2), false);
xlim([0 400])
ylim([0 90])
clim([0 7e-5])
title(['\DeltaOrientation pref.: ' exp{2} ' (n = ' num2str(n) ')'])
io.saveFigure(gcf, fPlots, sprintf('oriPref_%s', exp{2}))