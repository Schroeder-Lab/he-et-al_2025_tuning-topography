function Figure06(folders)

%% Parameters
minEV = 0.01; % minimum explained variance to plot RF
minPeak = 5; % minimum peak of RF (compared to noise) to plot RF

%% Examples

%% For all plots
fPlots = fullfile(folders.plots, 'Figures', 'Figure06');
if ~isfolder(fPlots)
    mkdir(fPlots)
end

%% Example spike waveshapes, firing traces and tuning curves

%% Tuning preferences and selectivity against depth in SC
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

        % select units with good RFs from noise or circles
        results = cell(2,1);
        if isfile(fullfile(f, '_ss_sparseNoise.times.npy'))
            results{1} = io.getNoiseRFFits(f);
        end
        if isfile(fullfile(f, 'circles.times.npy'))
            results{2} = io.getCircleRFFits(f);
        end
        stimTypes = rf.selectRFStim(results, minEV, minPeak);
        
        unitsSC = find(spikeData.clusterDepths(:,2) > 0);
        for k = 1:length(unitsSC)
            unit = unitsSC(k);
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


        chanCoord = readNPY(fullfile(f, 'channels.localCoordinates.npy'));
        SC_depth = readNPY(fullfile(f, '_ss_recordings.scChannels.npy'));
        SC_top = chanCoord(SC_depth(1), 2);
        SC_SO = chanCoord(SC_depth(2), 2);
    end
end

%% Population data: tuning preferences and selectivity against depth in SC

%% Pairwise differences in tuning preferences versus distance in depth

%% Example RFs

%% Compare preferred direction / orientation to expected value based on RF location