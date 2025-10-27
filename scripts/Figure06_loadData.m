function data = Figure06_loadData(folders, maxP, minEV, minPeak)

data = struct('animal', {}, 'date', {}, ...
    'dirPreferences', {}, 'oriPreferences', {}, ...
    'dirSel', {}, 'oriSel', {}, 'dirTuned', {}, 'oriTuned', {}, ...
    'depth', {}, 'totalN', {}, 'SO_depth', {}, ...
    'rfPos', {});

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
        rfData = cell(2,1);
        if isfile(fullfile(f, '_ss_sparseNoise.times.npy'))
            rfData{1} = io.getNoiseRFFits(f);
        end
        if isfile(fullfile(f, 'circles.times.npy'))
            rfData{2} = io.getCircleRFFits(f);
        end

        % select units with good RFs from noise or circles
        stimTypes = rf.selectRFStim(rfData, minEV, minPeak);
        
        unitsSC = spikeData.clusterDepths(:,2) > 0;
        if sum(unitsSC) < 1
            continue
        end
        data(session).animal = name;
        data(session).date = date;
        data(session).totalN = sum(unitsSC);
        SC_top = chanCoord(SC_depth(1), 2);
        SC_SO = chanCoord(SC_depth(2), 2);
        data(session).SO_depth = SC_top - SC_SO;

        responsive = unitsSC & ...
            (p_grandMean < maxP | p_compareAcrossStimuli < maxP);
        data(session).dirPreferences = ...
            dirTuning.preference(responsive);
        data(session).oriPreferences = ...
            oriTuning.preference(responsive);
        data(session).dirSel = dirTuning.selectivity(responsive);
        data(session).oriSel = oriTuning.selectivity(responsive);
        data(session).dirTuned = dirTuning.pValue(responsive) < maxP;
        data(session).oriTuned = oriTuning.pValue(responsive) < maxP;
        data(session).depth = spikeData.clusterDepths(responsive, 1);

        pos = NaN(sum(responsive), 2);
        for stim = 1:2
            rd = rfData{stim};
            if isempty(rd)
                continue
            end
            units_all = responsive & stimTypes == stim & ...
                rd.EV >= minEV & rd.peaks >= minPeak;
            units_tuning = units_all(responsive);
            pos(units_tuning,:) = rd.gaussPars(units_all, [2 4]);
        end
        data(session).rfPos = pos;

        session = session + 1;
    end
end

fprintf('    %d SC units were recorded across %d sessions from %d mice\n', ...
    sum([data.totalN]), length(data), ...
    length(unique({data.animal})))
numResponsive = length(cat(1,data.depth));
numDirTuned = sum(cat(1, data.dirTuned));
numOriTuned = sum(cat(1, data.oriTuned));
fprintf('    Of %d responsive units, %d (%.1f%%) were tuned to direction, %d (%.1f%%) were tuned to orientation\n', ...
    numResponsive, numDirTuned, ...
    numDirTuned / numResponsive * 100, ...
    numOriTuned, numOriTuned / numResponsive * 100)