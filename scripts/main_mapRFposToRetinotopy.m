% Find mapping between RF position and brain position (retinotopy). Then
% infer RF position of units where RF could not be mapped.

%% Folders
getFolders;

%% Parameters

%% Add paths
addpath(genpath(fullfile(folders.tools, 'npy-matlab')))
addpath(fullfile(folders.repo))

%% Fit RFs and get cross-validated explained variance
sets = {'boutons', 'neurons'};
for s = 1:2 % boutons and neurons
    fPlots = fullfile(folders.plots, 'RetinotopicRFs', sets{s});
    if ~isfolder(fPlots)
        mkdir(fPlots)
    end

    subjDirs = dir(fullfile(folders.data, sets{s}, 'SS*'));
    for subj = 1:length(subjDirs) % animals
        name = subjDirs(subj).name;
        fprintf('%s\n', name)
        dateDirs = dir(fullfile(folders.data, sets{s}, name, '2*'));
        for dt = 1:length(dateDirs) %dates
            date = dateDirs(dt).name;
            fprintf('  %s\n', date)
            f = fullfile(folders.data, sets{s}, name, date);
            % ignore session if visual noise stimuolus was not present
            if ~isfile(fullfile(f, '_ss_sparseNoise.times.npy'))
                continue
            end

            % load data
            data = io.getRecordingInfo(f);
            brainPos = data.roiPositions;
            data = io.getVisNoiseInfo(f);
            edges = data.edges;
            

            data = io.getCalciumData(f);
            time_tr = data.time;
            tr = data.traces;
            planes = data.planes;
            cellIDs = data.ids;
            delays = data.delays;
            data = io.getVisNoiseInfo(f);
            time_stim = data.times;
            stimMaps = data.frames;
            stimSeq = data.stimOrder;
            stimPos = data.edges; % [left right top bottom]
            % flip sign of stimulus borders along y-axis -> positive numbers
            % are above horizon/monitor centre
            stimPos(3:4) = -stimPos(3:4);
