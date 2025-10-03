function Figure06(folders)

%% Parameters

%% Examples

%% For all plots
fPlots = fullfile(folders.plots, 'Figures', 'Figure06');
if ~isfolder(fPlots)
    mkdir(fPlots)
end

%% Example spike waveshapes, firing traces and tuning curves

%% Example recordings: tuning preferences and selectivity against depth in SC
subjDirs = dir(fullfile(folders.data, 'ephys'));
subjDirs = subjDirs(~startsWith({subjDirs.name}, '.') & [subjDirs.isdir]);
for subj = 1:length(subjDirs) % animals
    name = subjDirs(subj).name;
    dateDirs = dir(fullfile(folders.data, 'ephys', name, '2*'));
    for dt = 1:length(dateDirs) %dates
        date = dateDirs(dt).name;
        f = fullfile(folders.data, 'ephys', name, date);

        % load data
        spikeData = io.getEphysData(f);
        stimData = io.getGratingInfo(f, 'gratingsDrifting');
        chanCoord = readNPY(fullfile(f, 'channels.localCoordinates.npy'));
        SC_depth = readNPY(fullfile(f, '_ss_recordings.scChannels.npy'));
        SC_top = chanCoord(SC_depth(1), 2);
        SC_SO = chanCoord(SC_depth(2), 2);


        % select units with good RFs from noise or circles
        [units, stimTypes] = rf.selectRFStim(results, minEV, minPeak);
        
        % determine depth of each unit
        depths = NaN(length(units), 1);
        for iUnit = 1: length(units)
            depths(iUnit) = median(spikeData.depths(...
                spikeData.clusters == units(iUnit)), "omitnan");
        end
        depths = SC_top - depths;
        % only consider units within SC
        ind = depths >= 0 & depths <= SC_extent;
        depths = depths(ind);
        units = units(ind);

        for k = 1:length(units)
            unit = units(k);
            res = results{stimTypes(k)};
            iUnit = find(res.units == unit);
            pars = res.fitParameters(iUnit,:);
            % ellipse at 2 STD (x and y), not rotated, not shifted
            x = pars(3) * cos(ellipse_x);
            y = pars(5) * sin(ellipse_x);
            % rotate and shift ellipse
            x_rot = pars(2) + x .* cos(pars(6)) - y .* sin(pars(6));
            y_rot = pars(4) + x .* sin(pars(6)) + y .* cos(pars(6));
            plot(x_rot, y_rot, lines{stimTypes(k)}, ...
                'Color', colors(colInds(k),:))
        end
    end
end

%% Population data: tuning preferences and selectivity against depth in SC

%% Pairwise differences in tuning preferences versus distance in depth

%% Example RFs

%% Compare preferred direction / orientation to expected value based on RF location