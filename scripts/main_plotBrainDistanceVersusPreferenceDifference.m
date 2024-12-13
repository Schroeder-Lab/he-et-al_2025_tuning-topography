function main_plotBrainDistanceVersusPreferenceDifference(folders)
% Plot pairwise differences in preferred direction or orientation against
% the distance of ROIs in the brain. Consider only ROIs within the same
% plane or ROIs across all planes (ignoring distances in depth).

%% Parameters
maxP = 0.05; % p-value threshold for response kernel and 
             % direction/orientation selectivity
minROIs = 15;
binSize = [5, 20];
stepSize = [2.5, 5];
numPerm = 1000;

%% Plot pairwise distance in brain versus difference in tuning preference
sets = {'boutons', 'neurons'};
for s = 1:2 % neurons and boutons
    subjDirs = dir(fullfile(folders.data, sets{s}, 'SS*'));
    dirDistPlane = {};
    dirDiffPlane = {};
    dirDiffPlaneNull = {};
    dirDistRec = {};
    dirDiffRec = {};
    dirDiffRecNull = {};

    oriDistPlane = {};
    oriDiffPlane = {};
    oriDiffPlaneNull = {};
    oriDistRec = {};
    oriDiffRec = {};
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

            % loop across each plane of recording + all planes together
            uniquePlanes = unique(planes);
            indPlanes = true(length(planes), length(uniquePlanes)+1);
            for p = 1:length(uniquePlanes)
                indPlanes(:,p) = planes == uniquePlanes(p);
            end
            dirDistPlane{rec} = [];
            dirDiffPlane{rec} = [];
            dirDiffPlaneNull{rec} = [];
            dirDistRec{rec} = [];
            dirDiffRec{rec} = [];
            dirDiffRecNull{rec} = [];

            oriDistPlane{rec} = [];
            oriDiffPlane{rec} = [];
            oriDiffPlaneNull{rec} = [];
            oriDistRec{rec} = [];
            oriDiffRec{rec} = [];
            oriDiffRecNull{rec} = [];
            for p = 1:size(indPlanes,2)
                % get tuning preferences of ROIs within plane
                validDir = indPlanes(:,p) & dirTuning.pValue < maxP;
                % for all unit pairs, determine distance in brain (ignore
                % depth);
                ddist = spatial.determineDistance(roiPos(validDir,1), ...
                    roiPos(validDir,2));
                dp = dirTuning.preference(validDir);
                % for all unit pairs, determine difference between preferred
                % directions/orientations
                ddiff = tuning.determinePreferenceDiff(dp, 'dir');
                % permute preferences to test significance
                ddiffPermuted = NaN(length(ddiff), numPerm);
                rng('default');
                for k = 1:numPerm
                    order = randperm(length(dp));
                    ddiffPermuted(:,k) = tuning.determinePreferenceDiff( ...
                        dp(order), 'dir');
                end

                validOri = indPlanes(:,p) & oriTuning.pValue < maxP;
                % for all unit pairs, determine distance in brain (ignore
                % depth);
                odist = spatial.determineDistance(roiPos(validOri,1), ...
                    roiPos(validOri,2));
                op = oriTuning.preference(validOri);
                % for all unit pairs, determine difference between preferred
                % orientations
                odiff = tuning.determinePreferenceDiff(op, 'ori');
                % permute preferences to test significance
                odiffPermuted = NaN(length(odiff), numPerm);
                rng('default');
                for k = 1:numPerm
                    order = randperm(length(op));
                    odiffPermuted(:,k) = tuning.determinePreferenceDiff( ...
                        op(order), 'ori');
                end

                % collect results
                if p < size(indPlanes,2) % for each plane
                    dirDistPlane{rec} = [dirDistPlane{rec}; ddist];
                    dirDiffPlane{rec} = [dirDiffPlane{rec}; ddiff];
                    dirDiffPlaneNull{rec} = [dirDiffPlaneNull{rec}; ddiffPermuted];
                    oriDistPlane{rec} = [oriDistPlane{rec}; odist];
                    oriDiffPlane{rec} = [oriDiffPlane{rec}; odiff];
                    oriDiffPlaneNull{rec} = [oriDiffPlaneNull{rec}; odiffPermuted];
                else % across all planes
                    dirDistRec{rec} = [dirDistRec{rec}; ddist];
                    dirDiffRec{rec} = [dirDiffRec{rec}; ddiff];
                    dirDiffRecNull{rec} = [dirDiffRecNull{rec}; ddiffPermuted];
                    oriDistRec{rec} = [oriDistRec{rec}; odist];
                    oriDiffRec{rec} = [oriDiffRec{rec}; odiff];
                    oriDiffRecNull{rec} = [oriDiffRecNull{rec}; odiffPermuted];
                end
                
                % plot results per recording (each + across planes)
                if p < size(indPlanes,2) % for each plane
                    strDir = sprintf('Direction_Plane%02d.jpg', uniquePlanes(p));
                    strOri = sprintf('Orientation_Plane%02d.jpg', uniquePlanes(p));
                else
                    strDir = 'Direction_allPlanes.jpg';
                    strOri = 'Orientation_allPlanes.jpg';
                end
                fig = spatial.plotPrefDiffVsDist(ddist, ddiff, ddiffPermuted, ...
                    binSize(s), stepSize(s), true);
                if ~isempty(fig)
                    title('\DeltaDirection pref. vs \Deltaposition')
                    saveas(gcf, fullfile(fPlots, strDir))
                    close gcf
                end
                fig = spatial.plotPrefDiffVsDist(odist, odiff, odiffPermuted, ...
                    binSize(s), stepSize(s), true);
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
    fig = spatial.plotPrefDiffVsDist(cat(1, dirDistPlane{:}), ...
        cat(1, dirDiffPlane{:}), cat(1, dirDiffPlaneNull{:}), ...
        binSize(s), stepSize(s), false);
    if ~isempty(fig)
        title('\DeltaDirection pref. vs \Deltaposition (within planes)')
        saveas(gcf, fullfile(fPlots, 'Direction_withinPlanes.jpg'))
        close gcf
    end
    fig = spatial.plotPrefDiffVsDist(cat(1, oriDistPlane{:}), ...
        cat(1, oriDiffPlane{:}), cat(1, oriDiffPlaneNull{:}), ...
        binSize(s), stepSize(s), false);
    if ~isempty(fig)
        title('\DeltaOrientation pref. vs \Deltaposition (within planes)')
        saveas(gcf, fullfile(fPlots, 'Orientation_withinPlanes.jpg'))
        close gcf
    end
    fig = spatial.plotPrefDiffVsDist(cat(1, dirDistRec{:}), ...
        cat(1, dirDiffRec{:}), cat(1, dirDiffRecNull{:}), ...
        binSize(s), stepSize(s), false);
    if ~isempty(fig)
        title('\DeltaDirection pref. vs \Deltaposition (across planes)')
        saveas(gcf, fullfile(fPlots, 'Direction_acrossPlanes.jpg'))
        close gcf
    end
    fig = spatial.plotPrefDiffVsDist(cat(1, oriDistRec{:}), ...
        cat(1, oriDiffRec{:}), cat(1, oriDiffRecNull{:}), ...
        binSize(s), stepSize(s), false);
    if ~isempty(fig)
        title('\DeltaOrientation pref. vs \Deltaposition (across planes)')
        saveas(gcf, fullfile(fPlots, 'Orientation_acrossPlanes.jpg'))
        close gcf
    end
end