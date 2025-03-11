function data = Figure04_loadData(folders, sets, retinotopyRF)

%% Parameters
% for evaluation of receptive fields (significance/goodness)
minEV = 0.01; % minimum explained variance to plot RF
minPeak = 5; % minimum peak of RF (compared to noise) to plot RF
dist2edge = 5; % minimum distance of RF centre to monitor edge
% for evaluation of direction/orientation selectivity
maxP = 0.05;

%% Load data: RF position, tuning preferences
data = struct('rfPos', cell(2,1), 'dirPref', [], 'DSI', [], ...
    'oriPref', [], 'OSI', [], 'set', []);
for s = 1:2
    subjDirs = dir(fullfile(folders.data, sets{s}, 'SS*'));
    count = 1;
    for subj = 1:length(subjDirs) % animals
        name = subjDirs(subj).name;
        dateDirs = dir(fullfile(folders.data, sets{s}, name, '2*'));
        for d = 1:length(dateDirs) %dates
            date = dateDirs(d).name;
            f = fullfile(folders.data, sets{s}, name, date);
            % ignore session if stimulus was not presented
            if ~isfile(fullfile(f, '_ss_gratingsDrifting.intervals.npy')) || ...
                    ~isfile(fullfile(f, '_ss_sparseNoise.times.npy'))
                continue
            end
            % load data
            if ~retinotopyRF(s)
                dt = io.getRFFits(f);
                edges = dt.edges; % [left right top bottom]
                pos = dt.fitParameters(:,[2 4]);
                % exclude all RFs too close to monitor edge
                ind = dt.EV < minEV | dt.peaks < minPeak | ...
                    pos(:,1) < edges(1)+dist2edge | ...
                    pos(:,1) > edges(2)-dist2edge | ...
                    pos(:,2) > edges(3)-dist2edge | ...
                    pos(:,2) < edges(4)+dist2edge;
                pos(ind,:) = NaN;
                clear dt
            else
                if ~isfile(fullfile(f, '_ss_rf.posRetinotopy.npy'))
                    continue
                end
                pos = readNPY(fullfile(f, '_ss_rf.posRetinotopy.npy'));
                edges = readNPY(fullfile(f, '_ss_rfDescr.edges.npy'));
                % exclude all RFs too close to monitor edge
                ind = pos(:,1) < edges(1)+dist2edge | ...
                    pos(:,1) > edges(2)-dist2edge | ...
                    pos(:,2) > edges(3)-dist2edge | ...
                    pos(:,2) < edges(4)+dist2edge;
                pos(ind,:) = NaN;
            end
            data(s).rfPos = [data(s).rfPos; pos];
            [dirTuning, oriTuning] = io.getTuningResults(f, 'gratingsDrifting');
            dp = dirTuning.preference;
            dsi = dirTuning.selectivity;
            dp(dirTuning.pValue >= maxP) = NaN;
            dsi(dirTuning.pValue >= maxP) = NaN;
            op = oriTuning.preference;
            osi = oriTuning.selectivity;
            op(oriTuning.pValue >= maxP) = NaN;
            osi(oriTuning.pValue >= maxP) = NaN;
            data(s).dirPref = [data(s).dirPref; dp];
            data(s).DSI = [data(s).DSI; dsi];
            data(s).oriPref = [data(s).oriPref; op];
            data(s).OSI = [data(s).OSI; osi];
            data(s).set = [data(s).set; ones(size(dp)) .* count];

            count = count + 1;
        end
    end
end