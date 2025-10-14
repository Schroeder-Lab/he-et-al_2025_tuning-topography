function Figure06(folders)

%% Parameters
minEV = 0.01; % minimum explained variance to plot RF
minPeak = 5; % minimum peak of RF (compared to noise) to plot RF
maxP = 0.05; % p-value threshold

%% Examples

%% For all plots
fPlots = fullfile(folders.plots, 'Figures', 'Figure06');
if ~isfolder(fPlots)
    mkdir(fPlots)
end

%% Example spike waveshapes, firing traces and tuning curves

%% Assemble tuning data and depth in SC across datasets
% variables of size #units
dirPreferences = [];
oriPreferences = [];
dirSel = [];
oriSel = [];
dirTuned = false(0,0);
oriTuned = false(0,0);
depth = []; % in um, relative to SC surface
dataset = [];
% variables of size #sessions
totalN = [];
animals = cell(0,1);
dates = cell(0,1);
SO_depth = [];

subjDirs = dir(fullfile(folders.data, 'ephys'));
subjDirs = subjDirs(~startsWith({subjDirs.name}, '.') & [subjDirs.isdir]);
session = 1;
for subj = 1:length(subjDirs) % animals
    name = subjDirs(subj).name;
    dateDirs = dir(fullfile(folders.data, 'ephys', name, '2*'));
    for dt = 1:length(dateDirs) %dates
        date = dateDirs(dt).name;
        f = fullfile(folders.data, 'ephys', name, date);
        % ignore session if stimulus was not presented
        if ~isfile(fullfile(f, '_ss_gratingsDrifting.intervals.npy'))
            continue
        end

        % load data
        spikeData = io.getEphysData(f);
        p_grandMean = readNPY(fullfile(f, ...
            '_ss_gratingsDriftingResponsive.p_grandMean.npy'));
        p_compareAcrossStimuli = readNPY(fullfile(f, ...
            '_ss_gratingsDriftingResponsive.p_acrossStimuli.npy'));
        [dirTuning, oriTuning] = io.getTuningResults(f, 'gratingsDrifting');
        chanCoord = readNPY(fullfile(f, 'channels.localCoordinates.npy'));
        SC_depth = readNPY(fullfile(f, '_ss_recordings.scChannels.npy'));
        
        unitsSC = spikeData.clusterDepths(:,2) > 0;
        if sum(unitsSC) < 1
            continue
        end
        animals{end+1,1} = name;
        dates{end+1,1} = date;
        totalN = [totalN; sum(unitsSC)];
        SC_top = chanCoord(SC_depth(1), 2);
        SC_SO = chanCoord(SC_depth(2), 2);
        SO_depth(end+1,1) = SC_top - SC_SO;

        responsive = unitsSC & ...
            (p_grandMean < maxP | p_compareAcrossStimuli < maxP);
        dirPreferences = [dirPreferences; dirTuning.preference(responsive)];
        oriPreferences = [oriPreferences; oriTuning.preference(responsive)];
        dirSel = [dirSel; dirTuning.selectivity(responsive)];
        oriSel = [oriSel; oriTuning.selectivity(responsive)];
        dirTuned = [dirTuned; dirTuning.pValue(responsive) < maxP];
        oriTuned = [oriTuned; oriTuning.pValue(responsive) < maxP];
        depth = [depth; spikeData.clusterDepths(responsive, 1)];
        dataset = [dataset; ones(sum(responsive),1) .* session];

        session = session + 1;
    end
end

fprintf('    %d SC units were recorded across %d sessions from %d mice\n', ...
    sum(totalN), length(totalN), length(unique(animals)))
fprintf('    Of %d responsive units, %d (%.1f%%) were tuned to direction, %d (%.1f%%) were tuned to orientation\n', ...
    length(depth), sum(dirTuned), sum(dirTuned)/length(depth)*100, ...
    sum(oriTuned), sum(oriTuned)/length(depth)*100)

%% Population data: tuning preferences and selectivity against depth in SC
scale = 30;
figure
hold on
for session = 1:length(totalN)
    units = find(dataset == session & dirTuned);
    X = ones(length(units),1) .* (session * 2 * scale);
    Y = depth(units);
    directions = dirPreferences(units);
    [U, V] = pol2cart(deg2rad(directions), ones(length(units),1) .* scale);
    X = X - 0.5 .* U;
    Y = Y - 0.5 .* V;
    quiver(X, Y, U, V, "off", 'k');
end

axis image
set(gca, "YDir", "reverse")


%% Pairwise differences in tuning preferences versus distance in depth

%% Example RFs

%% Compare preferred direction / orientation to expected value based on RF location


        % select units with good RFs from noise or circles
        results = cell(2,1);
        if isfile(fullfile(f, '_ss_sparseNoise.times.npy'))
            results{1} = io.getNoiseRFFits(f);
        end
        if isfile(fullfile(f, 'circles.times.npy'))
            results{2} = io.getCircleRFFits(f);
        end
        stimTypes = rf.selectRFStim(results, minEV, minPeak);