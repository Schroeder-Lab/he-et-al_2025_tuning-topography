% Plot the location of all units per recorded plane (using ROI masks) and
% color code preferred orientation or direction.

%% Folders
getFolders;

%% Parameters
maxP = 0.05;

types = {'dir','ori'};
strings = {'Direction','Orientation'};

%% Add paths
addpath(genpath(fullfile(folders.tools, 'npy-matlab')))
addpath(fullfile(folders.repo))

%% Plot direction and orientation maps
sets = {'neurons', 'boutons'};
for s = 1:2 % neurons and boutons
    subjDirs = dir(fullfile(folders.data, sets{s}, 'SS*'));
    for subj = 1:length(subjDirs) % animals
        name = subjDirs(subj).name;
        dateDirs = dir(fullfile(folders.data, sets{s}, name, '2*'));
        for dt = 1:length(dateDirs) %dates
            date = dateDirs(dt).name;
            f = fullfile(folders.data, sets{s}, name, date);
            % ignore session if stimulus was not presented
            if ~isfile(fullfile(f, '_ss_gratingsDrifting.intervals.npy'))
                continue
            end
            fPlots = fullfile(folders.plots, 'PreferenceMaps', ...
                sets{s}, name, date);
            if ~isfolder(fPlots)
                mkdir(fPlots)
            end

            % load data
            [dirTuning, oriTuning] = io.getTuningResults(f, 'gratingsDrifting');
            data = io.getCalciumData(f);
            planes = data.planes;
            data = io.getRecordingInfo(f);
            masks = data.roiMasks;
            fovPix = data.fovPix;
            fovM = data.fovMicrons;

            uniPlanes = unique(planes);
            for p = uniPlanes'
                indP = planes == p;
                tuning.plotOrientationMap(dirTuning.preference(indP), ...
                    dirTuning.pValue(indP) < maxP, 'dir', masks(indP,:), ...
                    fovPix(p,:), fovM(p,:));
                saveas(gcf, fullfile(fPlots, ...
                    sprintf('DirectionMap_Plane%02d.jpg', p)))
                close(gcf)
                tuning.plotOrientationMap(oriTuning.preference(indP), ...
                    oriTuning.pValue(indP) < maxP, 'ori', masks(indP,:), ...
                    fovPix(p,:), fovM(p,:));
                saveas(gcf, fullfile(fPlots, ...
                    sprintf('OrientationMap_Plane%02d.jpg', p)))
                close(gcf)
            end

        end
    end
end