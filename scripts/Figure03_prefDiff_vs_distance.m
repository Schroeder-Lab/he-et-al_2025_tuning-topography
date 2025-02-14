function Figure03_prefDiff_vs_distance(folders, sets, fPlots)

%% Parameters
% for evaluation of receptive fields (significance/goodness)
minEV = 0.01;
minPeak = 5;
% for evaluation of tuning selectivity
maxP = 0.05; % p-value threshold for direction/orientation selectivity
% for plotting
minROIs = 15;

% CONTINUE WITH: choose best binSize, stepSize, xLims for boutons and neurons
binSize = [1, 2];
stepSize = [0.2, 1];
xLims = [15 40];
% for testing
numPerm = 1000;

%% Plot pairwise distance in brain versus difference in tuning preference
for s = 1:2 % boutons and neurons
    subjDirs = dir(fullfile(folders.data, sets{s}, 'SS*'));
    dirDist = {};
    dirDiff = {};
    dirDiffNull = {};
    dirDiffRelative = {};
    rfDistBinnedDir = {};
    oriDist = {};
    oriDiff = {};
    oriDiffNull = {};
    oriDiffRelative = {};
    rfDistBinnedOri = {};
    rec = 1;
    for subj = 1:length(subjDirs) % animals
        name = subjDirs(subj).name;
        fprintf('%s\n', name)
        dateDirs = dir(fullfile(folders.data, sets{s}, name, '2*'));
        for dt = 1:length(dateDirs) %dates
            date = dateDirs(dt).name;
            f = fullfile(folders.data, sets{s}, name, date);
            % ignore session if stimulus was not presented
            if ~isfile(fullfile(f, '_ss_gratingsDrifting.intervals.npy')) || ...
                    ~isfile(fullfile(f, '_ss_rf.posRetinotopy.npy'))
                continue
            end
                
            % load data
            data = io.getRFFits(f);
            rfPos = data.fitParameters(:,[2 4]);
            ev_rf = data.EV;
            rf_peaks = data.peaks;
            [dirTuning, oriTuning] = io.getTuningResults(f, 'gratingsDrifting');

            dp = dirTuning.preference;
            valid = ~isnan(dp) & dirTuning.pValue < maxP & ...
                ev_rf >= minEV & rf_peaks >= minPeak;
            if sum(valid) < minROIs
                dirDist{rec} = [];
                dirDiff{rec} = [];
                dirDiffNull{rec} = [];
                dirDiffRelative{rec} = [];
                rfDistBinnedDir{rec} = [];
            else
                % for all unit pairs, determine distance of RFs
                ddist = spatial.determineDistance(rfPos(valid,1), ...
                    rfPos(valid,2));
                dp = dp(valid);
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
            valid = ~isnan(op) & oriTuning.pValue < maxP & ...
                ev_rf >= minEV & rf_peaks >= minPeak;
            if sum(valid) < minROIs
                oriDist{rec} = [];
                oriDiff{rec} = [];
                oriDiffNull{rec} = [];
                oriDiffRelative{rec} = [];
                rfDistBinnedOri{rec} = [];
            else
                % for all unit pairs, determine distance of RFs;
                odist = spatial.determineDistance(rfPos(valid,1), ...
                    rfPos(valid,2));
                op = op(valid);
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

    % plot across all datasets
    fig = spatial.plotPrefDiffVsDist(cat(1, dirDist{:}), ...
        cat(1, dirDiff{:}), cat(1, dirDiffNull{:}), ...
        binSize(s), stepSize(s), false);
    xlim([0 xLims(s)])
    ylim([0 180])
    xlabel('Distance (vis. deg.)')
    title('\DeltaDirection pref. vs \DeltaRF-position')
    io.saveFigure(fig, fPlots, ...
        sprintf('rfDistanceAll_%s_direction', sets{s}))
    fig = spatial.plotPrefDiffVsDist(cat(1, oriDist{:}), ...
        cat(1, oriDiff{:}), cat(1, oriDiffNull{:}), ...
        binSize(s), stepSize(s), false);
    xlim([0 xLims(s)])
    ylim([0 90])
    title('\DeltaOrientation pref. vs \DeltaRF-position')
    io.saveFigure(fig, fPlots, ...
        sprintf('rfDistanceAll_%s_orientation', sets{s}))
end