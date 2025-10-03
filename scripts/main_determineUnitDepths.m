function main_determineUnitDepths(folders)
% Given the site of the SC surface and its SGS-SO border, determine depth
% of each unit within SC.

%% Parameters
% for depth estimate
SC_extent = 1000; % microns

%% Determine depth of each unit
subjDirs = dir(fullfile(folders.data, 'ephys'));
subjDirs = subjDirs(~startsWith({subjDirs.name}, '.') & [subjDirs.isdir]);
for subj = 1:length(subjDirs) % animals
    name = subjDirs(subj).name;
    fprintf('%s\n', name)
    dateDirs = dir(fullfile(folders.data, 'ephys', name, '2*'));
    for dt = 1:length(dateDirs) %dates
        date = dateDirs(dt).name;
        fprintf('  %s\n', date)
        f = fullfile(folders.data, 'ephys', name, date);

        % Load data
        spikeData = io.getEphysData(f);
        chanCoord = readNPY(fullfile(f, 'channels.localCoordinates.npy'));
        SC_depth = readNPY(fullfile(f, '_ss_recordings.scChannels.npy'));
        SC_top = chanCoord(SC_depth(1), 2);
        SC_SO = chanCoord(SC_depth(2), 2);
        if isfile(fullfile(f, '_ss_sparseNoise.times.npy'))
            stimData = io.getVisNoiseInfo(f);
        elseif isfile(fullfile(f, 'circles.times.npy'))
            stimData = io.getCircleInfo(f);
        else
            fprintf('  NO NOISE OR CIRCLE DATA. DISREGARD DATASET!\n')
            continue
        end

        % determine depth of all units along probe, estimated across
        % complete recording
        depths = NaN(length(spikeData.clusterIDs), 2);
        % spike times within stimulus paradigm used to determine SC depth
        t_ind = spikeData.times >= stimData.times(1) - 5 & ...
            spikeData.times <= stimData.times(end) + 1;
        for k = 1:length(spikeData.clusterIDs)
            ind_unit = spikeData.clusters == spikeData.clusterIDs(k);
            if any(t_ind & ind_unit)
                depths(k,1) = median( ...
                    spikeData.depths(t_ind & ind_unit), "omitnan");
            else
                depths(k,1) = median(spikeData.depths(ind_unit), ...
                    "omitnan");
            end
        end

        % determine depth relative to SC surface
        depths(:,1) = SC_top - depths(:,1);
        % classify depths: (-1) above SC, (0) below SC, (1) sSC, (2) dSC
        depths(depths(:,1) < 0, 2) = -1;
        depths(depths(:,1) > SC_extent, 2) = 0;
        depths(depths(:,1) >= 0 & depths(:,1) <= SC_top - SC_SO, 2) = 1;
        depths(depths(:,1) > SC_top - SC_SO & ...
            depths(:,1) < SC_extent, 2) = 2;

        % save depth information
        writeNPY(depths, fullfile(f, 'clusters.depths.npy'))
    end
end