function Figure02SA(folders)

%% Parameters
sets = {'boutons', 'neurons'};
% for evaluation of receptive fields (significance/goodness)
minEV = 0.01;
minPeak = 5;
% to relate RF positions mapped with and without regards to eye movements
minUnits_RFs = 10;
% plotting RFs
edges_rf = [-135 -80 42 -20];
RFtypes = {'ON', 'OFF', 'ON+OFF'};

%% Examples
ex = cell(2,3); % rows: (1) bouton, (2) neuron
ex(1,:) = {'SS069', '2016-10-21', 12};
ex(2,:) = {'SS044', '2015-05-29', 7};

%% For all plots
fPlots = fullfile(folders.plots, 'Figures', 'Figure02SA');
if ~isfolder(fPlots)
    mkdir(fPlots)
end

%% Plot example eye data and RFs (corrected and not)
for s = 1:2
    name = ex{s,1};
    date = ex{s,2};
    f = fullfile(folders.data, sets{s}, name, date);

    data = io.getNoiseRFFits(f, true);

    % Eye positions
    A = data.A;
    gazePos = A * data.eyePos';

    figure
    subplot(2,1,1)
    scatter(gazePos(1,:), gazePos(2,:), 15, 'k', 'filled', 'o')
    axis image
    xlabel('\DeltaAzimuth (vis. deg.)')
    ylabel('\DeltaElevation (vis. deg.)')
    title(sprintf('%s %s', name, date))

    subplot(2,1,2)
    scatter(data.eyePos(:,1), data.eyePos(:,2), 15, 'k', 'filled', 'o')
    axis image
    set(gca, 'XDir', 'reverse')
    xlabel('\DeltaX (pixels)')
    ylabel('\DeltaY (pixels)')

    io.saveFigure(gcf, fPlots, ...
        sprintf('eye-positions_%s_scatter', sets{s}))

    % example RF
    figs = [0 0];
    data_uncorrected = io.getNoiseRFFits(f, false);
    edges = data.edges;
    gridW = diff(edges(1:2)) / size(data.maps,3);
    gridH = -diff(edges(3:4)) / size(data.maps,2);
    iUnit = ex{s,3};
    % uncorrected RF
    % rf_tmp: [rows x columns x subfield]
    rf_tmp = squeeze(mean(data_uncorrected.maps(iUnit,:,:,:,:),4));
    rf.plotRF(rf_tmp, data_uncorrected.gaussPars(iUnit,:), ...
        data_uncorrected.bestSubFields(iUnit), edges, edges_rf, gridW, gridH)
    sgtitle(sprintf('ROI %d (EV: %.3f, peak/noise: %.1f, %s)', ...
        iUnit, data_uncorrected.EV(iUnit), data_uncorrected.peaks(iUnit), ...
        RFtypes{data_uncorrected.bestSubFields(iUnit)}))
    lim1 = clim;
    figs(1) = gcf;

    % corrected RF
    % rf_tmp: [rows x columns x subfield]
    rf_tmp = squeeze(mean(data.maps(iUnit,:,:,:,:),4));
    rf.plotRF(rf_tmp, data.gaussPars(iUnit,:), ...
        data.bestSubFields(iUnit), edges, edges_rf, gridW, gridH)
    sgtitle(sprintf('ROI %d (EV: %.3f, peak/noise: %.1f, %s)', ...
        iUnit, data.EV(iUnit), data.peaks(iUnit), ...
        RFtypes{data.bestSubFields(iUnit)}))
    lim2 = clim;
    figs(2) = gcf;

    limit_x = [-1 1] .* max(lim1(2), lim2(2));
    for k = 1:2
        figure(figs(k))
        for j = 1:2
            subplot(1,2,j)
            clim(limit_x)
        end
    end

    figure(figs(1))
    io.saveFigure(gcf, fPlots, sprintf('example_%s_RF-uncorrected_%s_%s_%03d', ...
        sets{s}, ex{s,1}, ex{s,2}, iUnit))
    figure(figs(2))
    io.saveFigure(gcf, fPlots, sprintf('example_%s_RF-corrected_%s_%s_%03d', ...
        sets{s}, ex{s,1}, ex{s,2}, iUnit))
end

%% Plot distance between RF centres mapped with vs without regards to eye movements
for s = 1:2 % boutons and neurons
    subjDirs = dir(fullfile(folders.data, sets{s}, 'SS*'));
    rfPos = [];
    rfPos_eye = [];

    % Collect data
    for subj = 1:length(subjDirs) % animals
        name = subjDirs(subj).name;
        dateDirs = dir(fullfile(folders.data, sets{s}, name, '2*'));
        for dt = 1:length(dateDirs) %dates
            date = dateDirs(dt).name;
            f = fullfile(folders.data, sets{s}, name, date);
            % ignore session if visual noise stimuolus was not present
            if ~isfile(fullfile(f, '_ss_rf_eye.gaussFitPars.npy'))
                continue
            end

            data = io.getNoiseRFFits(f, false);
            valid = data.EV >= minEV & data.peaks >= minPeak;
            pos = data.gaussPars(:, [2, 4]);

            data = io.getNoiseRFFits(f, true);
            valid_eye = data.EV >= minEV & data.peaks >= minPeak;
            pos_eye = data.gaussPars(:, [2, 4]);

            valid = valid & valid_eye;
            if sum(valid) < minUnits_RFs
                continue
            end

            rfPos = [rfPos; pos(valid,:)];
            rfPos_eye = [rfPos_eye; pos_eye(valid,:)];
        end
    end

    % Distance of eye-corrected RFs from not-corrected RFs
    dist = rfPos_eye - rfPos;

    euclidDist = sqrt(sum(dist .^ 2, 2));
    euclidDist(euclidDist > 4) = 4;

    figure
    % density of data points
    dens = ksdensity(dist, dist);
    % scatterplot
    scatter(dist(:,1), dist(:,2), 10, dens, "filled")
    colormap(gca, flip(gray))
    limit_x = clim;
    clim([0 limit_x(2)])
    c = colorbar;
    c.Label.String = 'Density';
    axis image
    xlim([-1 1] .* 2)
    ylim([-1 1] .* 2)
    set(gca, 'box', 'off')
    xlabel('\DeltaAzimuth (vis. deg.)')
    ylabel('\DeltaElevation (vis.deg.)')
    title(sprintf('%s (n = %d)', sets{s}, length(dist)))
    io.saveFigure(gcf, fPlots, ...
        sprintf('rf-eye-corrected_%s_scatter', sets{s}))

    figure
    histogram(euclidDist, 0:0.2:4)
    set(gca, 'box', 'off')
    xlabel('Euclidean distance (vis. deg.)')
    ylabel(sprintf('#%s', sets{s}))
    title(sprintf('Median: %.2f (%d%% < 1 deg, %d%% < 2 deg)', ...
        median(euclidDist, "omitnan"), ...
        round(sum(euclidDist < 1) / sum(~isnan(euclidDist)) * 100), ...
        round(sum(euclidDist < 2) / sum(~isnan(euclidDist)) * 100)))
    io.saveFigure(gcf, fPlots, ...
        sprintf('rf-eye-corrected_%s_histogram', sets{s}))
end

%% Plot RF position during gratings relative to the position mapped with noise
% one scatter per direction, pooled across bouton and neuron data
limit_x = 8;
limit_y = 4;
eyePos = cell(12,1); % eye position per grating direction
nSessions = 0;
for s = 1:2 % boutons and neurons
    subjDirs = dir(fullfile(folders.data, sets{s}, 'SS*'));
    for subj = 1:length(subjDirs) % animals
        name = subjDirs(subj).name;
        dateDirs = dir(fullfile(folders.data, sets{s}, name, '2*'));
        for dt = 1:length(dateDirs) %dates
            date = dateDirs(dt).name;
            f = fullfile(folders.data, sets{s}, name, date);
            % ignore session if eye data not present
            if ~isfile(fullfile(f, 'eyeToGaze.transformMatrix.npy')) || ...
                    ~isfile(fullfile(f, '_ss_gratingsDrifting.intervals.npy'))
                continue
            end

            rfData = io.getNoiseRFFits(f, true);
            A = rfData.A;
            defaultPos = rfData.defaultPos;
            eyeData = io.getPupilData(f);
            stimData = io.getGratingInfo(f, 'gratingsDrifting');

            valid = rfData.EV >= minEV & rfData.peaks >= minPeak;
            if sum(valid) < minUnits_RFs
                continue
            end

            % get eye position relative to default position (during visual
            % noise)
            pos = eyeData.xyPos - defaultPos';

            % transform eye position from pixels to visual degrees
            pos = (A * pos')';

            % select eye positions during grating presentations
            stimDur = min(diff(stimData.times, 1, 2));
            % pos: [t x trials x 2]
            pos = traces.getAlignedTraces(pos, eyeData.time, ...
                stimData.times(:,1), [0 stimDur]);

            for stim = 1:12
                ind = stimData.ids == stim;
                pos_stim = pos(:,ind,:);
                eyePos{stim} = [eyePos{stim}; reshape(pos_stim, [], 2)];
            end

            nSessions = nSessions + 1;
        end
    end
end

% Kernel density estimate on a grid
x = linspace(-limit_x, limit_x, 200);
y = linspace(-limit_y, limit_y, 100);
[X,Y] = meshgrid(x,y);
gridPts = [X(:), Y(:)];

plots = cell(12,1);
limit_z = 0;
for stim = 1:12
    % Evaluate 2-D ksdensity -> pdf values at grid points
    f = ksdensity(eyePos{stim}, gridPts, 'Function', 'pdf');  % returns f as column vector

    % Reshape and plot heatmap
    plots{stim} = reshape(f, size(X));
    m = max(plots{stim}, [], "all");
    if m > limit_z
        limit_z = m;
    end
end

for stim = 1:12
    plots{stim} = plots{stim} ./ limit_z;
    plots{stim} = 1 - plots{stim};
    imwrite(plots{stim}, fullfile(fPlots, ...
        sprintf('eyePosPerGrating_%d.png', stimData.directions(stim))));
end

fprintf("%d sessions used for eye movements during drifting gratings.\n", nSessions)