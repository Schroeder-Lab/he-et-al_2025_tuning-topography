%% Folders
getFolders;

%% Parameters
binSize = [10, 5];
stepSize = [5, 2.5];
numPerm = 1000;

%% Examples
ex = cell(2,3); % rows: (1) bouton, (2) neuron
ex(1,:) = {'SS078', '2017-09-28', 1};
ex(2,:) = {'SS044', '2015-04-28', 3};

%% Plot pairwise distance in brain versus difference in tuning preference
sets = {'neurons', 'boutons'};
stimTypes = {'gratingsDrifting', 'gratingsStatic', 'bars'};
for s = 2 % neurons and boutons
    subjDirs = dir(fullfile(folders.data, sets{s}, 'SS*'));
    distPlane = {};
    distRec = {};
    dirDiffPlane = {};
    dirDiffRec = {};
    dirDiffPlaneNull = {};
    dirDiffRecNull = {};
    oriDiffPlane = {};
    oriDiffRec = {};
    oriDiffPlaneNull = {};
    oriDiffRecNull = {};
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
            fPlots = fullfile(folders.plots, ...
                'PreferenceDiffVsBrainDistance', sets{s}, name, date);
            if ~isfolder(fPlots)
                mkdir(fPlots)
            end
                
            % load data
            data = io.getCalciumData(f);
            planes = data.planes;
            data = io.getRecordingInfo(f);
            roiPos = data.roiPositions(:,1:2);
            [dirTuning, oriTuning] = io.getTuningResults(f, 'gratingsDrifting');

            uniquePlanes = unique(planes);
            indPlanes = true(length(planes), length(uniquePlanes)+1);
            for p = 1:length(uniquePlanes)
                indPlanes(:,p) = planes == uniquePlanes(p);
            end
            distPlane{rec} = [];
            dirDiffPlane{rec} = [];
            oriDiffPlane{rec} = [];
            distRec{rec} = [];
            dirDiffRec{rec} = [];
            oriDiffRec{rec} = [];
            dirDiffPlaneNull{rec} = [];
            dirDiffRecNull{rec} = [];
            oriDiffPlaneNull{rec} = [];
            oriDiffRecNull{rec} = [];
            for p = 1:size(indPlanes,2)
                dp = dirTuning.preference(indPlanes(:,p));
                op = oriTuning.preference(indPlanes(:,p));
                indValid1 = ~(isnan(dp) & isnan(op));
                dp = dp(indValid1);
                op = op(indValid1);
                indValid2 = find(indPlanes(:,p));
                indValid2 = indValid2(indValid1);

                % for all unit pairs, determine difference between preferred
                % directions/orientations
                dd = tuning.determinePreferenceDiff(dp, 'dir');
                od = tuning.determinePreferenceDiff(op, 'ori');
                % for all unit pairs, determine distance in brain (ignore
                % depth);
                d = spatial.determineDistance(roiPos(indValid2,1), ...
                    roiPos(indValid2,2));
                % permute preferences to test significance
                ddPermuted = NaN(length(dd), numPerm);
                odPermuted = NaN(length(od), numPerm);
                rng('default');
                for k = 1:numPerm
                    order = randperm(length(dp));
                    ddPermuted(:,k) = tuning.determinePreferenceDiff(dp(order), 'dir');
                    odPermuted(:,k) = tuning.determinePreferenceDiff(op(order), 'ori');
                end
                % collect results
                if p < size(indPlanes,2) % for each plane
                    distPlane{rec} = [distPlane{rec}; d];
                    dirDiffPlane{rec} = [dirDiffPlane{rec}; dd];
                    oriDiffPlane{rec} = [oriDiffPlane{rec}; od];
                    dirDiffPlaneNull{rec} = [dirDiffPlaneNull{rec}; ddPermuted];
                    oriDiffPlaneNull{rec} = [oriDiffPlaneNull{rec}; odPermuted];
                else % across all planes
                    distRec{rec} = [distRec{rec}; d];
                    dirDiffRec{rec} = [dirDiffRec{rec}; dd];
                    oriDiffRec{rec} = [oriDiffRec{rec}; od];
                    dirDiffRecNull{rec} = [dirDiffRecNull{rec}; ddPermuted];
                    oriDiffRecNull{rec} = [oriDiffRecNull{rec}; odPermuted];
                end
                
                % plot
                if p < size(indPlanes,2) % for each plane
                    strDir = sprintf('Direction_Plane%02d.jpg', uniquePlanes(p));
                    strOri = sprintf('Orientation_Plane%02d.jpg', uniquePlanes(p));
                else
                    strDir = 'Direction_allPlanes.jpg';
                    strOri = 'Orientation_allPlanes.jpg';
                end
                fig = spatial.plotPrefDiffVsDist(d, dd, ddPermuted, ...
                    binSize(s), stepSize(s));
                if ~isempty(fig)
                    title('\DeltaDirection pref. vs \Deltaposition')
                    saveas(gcf, fullfile(fPlots, strDir))
                    close gcf
                end
                fig = spatial.plotPrefDiffVsDist(d, od, odPermuted, ...
                    binSize(s), stepSize(s));
                if ~isempty(fig)
                    title('\DeltaOrientation pref. vs \Deltaposition')
                    saveas(gcf, fullfile(fPlots, strOri))
                    close gcf
                end
            end
            rec = rec + 1;
        end
    end

    % plot across all datasets
    fPlots = fullfile(folders.plots, ...
        'PreferenceDiffVsBrainDistance', sets{s});
    fig = spatial.plotPrefDiffVsDist(cat(1, distPlane{:}), ...
        cat(1, dirDiffPlane{:}), cat(1, dirDiffPlaneNull{:}), ...
        binSize(s), stepSize(s));
    if ~isempty(fig)
        title('\DeltaDirection pref. vs \Deltaposition (within planes)')
        saveas(gcf, fullfile(fPlots, 'Direction_withinPlanes.jpg'))
        close gcf
    end
    fig = spatial.plotPrefDiffVsDist(cat(1, distPlane{:}), ...
        cat(1, oriDiffPlane{:}), cat(1, oriDiffPlaneNull{:}), ...
        binSize(s), stepSize(s));
    if ~isempty(fig)
        title('\DeltaOrientation pref. vs \Deltaposition (within planes)')
        saveas(gcf, fullfile(fPlots, 'Orientation_withinPlanes.jpg'))
        close gcf
    end
    fig = spatial.plotPrefDiffVsDist(cat(1, distRec{:}), ...
        cat(1, dirDiffRec{:}), cat(1, dirDiffRecNull{:}), ...
        binSize(s), stepSize(s));
    if ~isempty(fig)
        title('\DeltaDirection pref. vs \Deltaposition (across planes)')
        saveas(gcf, fullfile(fPlots, 'Direction_acrossPlanes.jpg'))
        close gcf
    end
    fig = spatial.plotPrefDiffVsDist(cat(1, distRec{:}), ...
        cat(1, oriDiffRec{:}), cat(1, oriDiffRecNull{:}), ...
        binSize(s), stepSize(s));
    if ~isempty(fig)
        title('\DeltaOrientation pref. vs \Deltaposition (across planes)')
        saveas(gcf, fullfile(fPlots, 'Orientation_acrossPlanes.jpg'))
        close gcf
    end
end