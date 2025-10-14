function main_fitReceptiveFields_ephys(folders)
% Given the spike responses to the visual noise or circles stimulus, fit spatial
% receptive fields (ON and OFF fields).

%% Parameters
% for receptive field fits
lambda = 0.002; % smoothing parameter for spatial RF
rf_timeLimits = [0 0.2]; % time range of stimulus before neural response 
                      % considered for RF
RFtypes = {'ON', 'OFF', 'ON+OFF'};

% for evaluation of receptive fields (significance/goodness)
minEV = 0.01; % minimum explained variance to plot RF
minPeak = 5; % minimum peak of RF (compared to noise) to plot RF
minEV_plot = 0.005; % minimum explained variance to plot RF
minPeak_plot = 3.5; % minimum peak of RF (compared to noise) to plot RF

% plotting RFs
plotFolder = '05_ReceptiveFields';
% colormaps
[cm_ON, cm_OFF] = colmaps.getRFMaps;
cms = cat(3, cm_ON, cm_OFF);
ellipse_x = linspace(-pi, pi, 100);
titles = {'ON field','OFF field'};

% for depth estimate
SC_extent = 1000; % microns

%% Fit RFs and get cross-validated explained variance
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

        %% Load data
        spikeData = io.getEphysData(f);

        for stimType = 1:2
            if stimType == 1
                % if white noise data available, use it to map RFs
                if ~isfile(fullfile(f, '_ss_sparseNoise.times.npy'))
                    continue
                end
                stimData = io.getVisNoiseInfo(f);
                % edges: [left right top bottom] (above horizon: >0)
                % gridX, gridY: centre coordinates of checkers
                [stimMatrix, edges, gridX, gridY] = ...
                    stimuli.getNoiseStimMatrix(stimData.edges, ...
                    stimData.frames, stimData.stimOrder);
                stimSize = size(stimMatrix, [2 3]);
            else
                % if circle data available, use that to map RFs
                if ~isfile(fullfile(f, 'circles.times.npy'))
                    continue
                end
                stimData = io.getCircleInfo(f);
                [stimMatrix, gridX, gridY, diameters] = ...
                    circles.getStimMatrix(stimData.times, ...
                    stimData.xPos, stimData.yPos, stimData.diameter, ...
                    stimData.isWhite);
                stimSize = [length(gridY) length(gridX) length(diameters)];
                gridW = median(diff(gridX));
                gridH = median(-diff(gridY));
                % edges: [left right top bottom] (above horizon: >0)
                edges = [gridX(1)-0.5*gridW gridX(end)+0.5*gridW ...
                    gridY(1)+0.5*gridH gridY(end)-0.5*gridH];
            end
            t_stim = stimData.times;
            tBin_stim = median(diff(t_stim));
            rfBins = floor(rf_timeLimits(1) / tBin_stim) : ...
                ceil(rf_timeLimits(2) / tBin_stim) - 1;

            %% Prepare traces

            % ignore spike data before/after visual noise stimulation
            t_ind = spikeData.times >= stimData.times(1) - 10 & ...
                spikeData.times <= stimData.times(end) + 10;
            spikeData.times = spikeData.times(t_ind);
            spikeData.clusters = spikeData.clusters(t_ind);

            % get firing rates aligned to stimulus times
            traces = nan(length(t_stim), length(spikeData.clusterIDs));
            for iUnit = 1:length(spikeData.clusterIDs)
                t = spikeData.times(spikeData.clusters == ...
                    spikeData.clusterIDs(iUnit));
                if isempty(t)
                    continue
                end
                [~, frameOfSpike] = events.alignData(t, ...
                    t_stim, [0 tBin_stim]);
                traces(:,iUnit) = histcounts(frameOfSpike, ...
                    (0:length(t_stim)) + 0.5);
            end
            % z-score neural response
            zTraces = (traces - mean(traces,1,'omitnan')) ./ ...
                std(traces,0,1,'omitnan');

            %% Map RFs
            % generate toplitz matrix for stimulus
            [toeplitz, t_toeplitz] = ...
                rf.makeStimToeplitz(stimMatrix, t_stim, rfBins);

            %--------------------------------------------------------------
            % Comment if RFs are already mapped and only Gaussian fit is
            % needed.
            % map RF
            % rFields: [rows x columns (x diameters) x time x ON/OFF x units]
            rFields = rf.getReceptiveField( ...
                zTraces, toeplitz, stimSize, rfBins, lambda);
            % NOTE: OFF subfield: positive pixel -> suppressed by black
            % negative pixel -> driven by black

            %--------------------------------------------------------------
            % % Uncomment if RFs have already been mapped, and only Gaussian
            % % fit should be performed.
            % results = io.getRFFits(f);
            % rFields = results.maps;
            % rFields = permute(rFields, [2:5 1]);
            %--------------------------------------------------------------

            % fit Gaussian
            [rfGaussPars, fitGaussians, fitWeights, peakNoiseRatio, ...
                bestSubFields, subFieldSigns, predictions, EVs, sizeTuning] = ...
                rf.fitAllRFs(rFields, rfBins, gridX, gridY, zTraces, toeplitz);
            % save results
            dims = 1:ndims(rFields);
            results.maps = permute(rFields, dims([end 1:end-1]));
            results.bestSubfields = bestSubFields;
            results.subfieldSigns = subFieldSigns;
            results.gaussPars = rfGaussPars;
            results.peaks = peakNoiseRatio;
            results.gaussMasks = fitGaussians;
            results.timeWeights = fitWeights;
            results.EV = EVs;
            results.predictions = predictions;
            results.time_prediction = t_toeplitz;
            results.time_RF = rfBins * tBin_stim;
            if stimType == 1
                results.edges = edges;
                io.writeNoiseRFResults(results, f);
            elseif stimType == 2
                results.x = gridX;
                results.y = gridY;
                results.diameters = diameters;
                results.sizeTuning = sizeTuning;
                io.writeCircleRFResults(results, f);
            end
            clear results
        end
    end
end

%% Plot RFs
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

        spikeData = io.getEphysData(f);

        for stimType = 1:2
            if stimType == 1
                if ~isfile(fullfile(f, '_ss_sparseNoise.times.npy'))
                    continue
                end
                fPlots = fullfile(folders.plots, plotFolder, ...
                    'ephys', 'noise', name, date);
                if ~isfolder(fPlots)
                    mkdir(fullfile(fPlots, '1_above'))
                    mkdir(fullfile(fPlots, '2_sSC'))
                    mkdir(fullfile(fPlots, '3_dSC'))
                    mkdir(fullfile(fPlots, '4_below'))
                end
                % load data
                results = io.getNoiseRFFits(f);
                edges = results.edges;
                gridW = diff(edges(1:2)) / size(results.maps,3);
                gridH = -diff(edges(3:4)) / size(results.maps,2);
            else
                if ~isfile(fullfile(f, 'circles.times.npy'))
                    continue
                end
                fPlots = fullfile(folders.plots, plotFolder, ...
                    'ephys', 'circles', name, date);
                if ~isfolder(fPlots)
                    mkdir(fullfile(fPlots, '1_above'))
                    mkdir(fullfile(fPlots, '2_sSC'))
                    mkdir(fullfile(fPlots, '3_dSC'))
                    mkdir(fullfile(fPlots, '4_below'))
                end
                % load data
                results = io.getCircleRFFits(f);
                gridW = median(diff(results.x));
                gridH = median(-diff(results.y));
                % edges: [left right top bottom] (above horizon: >0)
                edges = [results.x(1)-0.5*gridW results.x(end)+0.5*gridW ...
                    results.y(1)+0.5*gridH results.y(end)-0.5*gridH];
            end

            for iUnit = 1:length(spikeData.clusterIDs)
                % plot RF (if explained var. and peak-to-noise high enough
                if isnan(results.EV(iUnit)) || ...
                        results.EV(iUnit) < minEV_plot || ...
                        results.peaks(iUnit) < minPeak_plot
                    continue
                end
                line = ':';
                if results.EV(iUnit) >= minEV && ...
                        results.peaks(iUnit) >= minPeak
                    line = '-';
                end
                [~, mxTime] = max(results.timeWeights(iUnit,:));
                % rf: [rows x columns x subfield]
                if ndims(results.maps) == 5 % visual noise
                    rfield = squeeze(results.maps(iUnit,:,:,mxTime,:));
                elseif ndims(results.maps) == 6 % circle paradigm
                    [~, mxSize] = max(results.sizeTuning(iUnit,:));
                    rfield = squeeze(results.maps(iUnit,:,:,mxSize,mxTime,:));
                end
                rfield(:,:,2) = -rfield(:,:,2);
                mx = max(abs(rfield),[],"all");
                rfGaussPars = results.gaussPars(iUnit,:);
                figure('Position', [75 195 1470 475])
                for sf = 1:2
                    subplot(1,2,sf)
                    % STA
                    imagesc([edges(1)+gridW/2 edges(2)-gridW/2], ...
                        [edges(3)-gridH/2 edges(4)+gridH/2], ...
                        rfield(:,:,sf),[-mx mx])
                    hold on
                    if ismember(results.bestSubFields(iUnit), [sf 3])
                        % ellipse at 2 STD (x and y), not rotated, not shifted
                        x = rfGaussPars(3) * cos(ellipse_x);
                        y = rfGaussPars(5) * sin(ellipse_x);
                        % rotate and shift ellipse
                        x_rot = rfGaussPars(2) + ...
                            x .* cos(rfGaussPars(6)) - ...
                            y .* sin(rfGaussPars(6));
                        y_rot = rfGaussPars(4) + ...
                            x .* sin(rfGaussPars(6)) + ...
                            y .* cos(rfGaussPars(6));
                        n = x_rot < edges(1) | x_rot > edges(2) | ...
                            y_rot > edges(3) | y_rot < edges(4);
                        x_rot(n) = NaN;
                        y_rot(n) = NaN;
                        plot(x_rot, y_rot, ['k' line], 'LineWidth', 2)
                    end
                    axis image
                    set(gca, 'box', 'off', 'YDir', 'normal')
                    xlim(edges([1 2]))
                    ylim(edges([4 3]))
                    colormap(gca, cms(:,:,sf))
                    title(titles{sf})
                    colorbar
                end

                depth = spikeData.clusterDepths(iUnit,:);
                if stimType == 1
                    sgtitle(sprintf(['ROI %d ' ...
                        '(EV: %.3f, peak/noise: %.1f, %s)\n' ...
                        '%.1f um from SC surface'], ...
                        spikeData.clusterIDs(iUnit), ...
                        results.EV(iUnit), results.peaks(iUnit), ...
                        RFtypes{results.bestSubFields(iUnit)}, ...
                        depth(1)))
                else
                    sgtitle(sprintf(['ROI %d, size: %.1f deg ' ...
                        '(EV: %.3f, peak/noise: %.1f, %s)\n' ...
                        '%.1f um from SC surface'], ...
                        spikeData.clusterIDs(iUnit), ...
                        results.diameters(mxSize), ...
                        results.EV(iUnit), results.peaks(iUnit), ...
                        RFtypes{results.bestSubFields(iUnit)}, ...
                        depth(1)))
                end

                switch depth(2)
                    case -1
                        fp = '1_above';
                    case 1
                        fp = '2_sSC';
                    case 2
                        fp = '3_dSC';
                    case 0
                        fp = '4_below';
                end
                saveas(gcf, fullfile(fPlots, fp, ...
                    sprintf('%04d_Unit%04d.jpg', round(depth(1)), ...
                    spikeData.clusterIDs(iUnit))));
                close gcf
            end
        end
    end
end

%% Plot all RF outlines per dataset
fPlots = fullfile(folders.plots, plotFolder, 'ephys');

subjDirs = dir(fullfile(folders.data, 'ephys'));
subjDirs = subjDirs(~startsWith({subjDirs.name}, '.') & [subjDirs.isdir]);
for subj = 1:length(subjDirs) % animals
    name = subjDirs(subj).name;
    dateDirs = dir(fullfile(folders.data, 'ephys', name, '2*'));
    for dt = 1:length(dateDirs) %dates
        date = dateDirs(dt).name;
        f = fullfile(folders.data, 'ephys', name, date);
        % ignore session if no visual noise or circle stimulus was present
        if ~isfile(fullfile(f, '_ss_sparseNoise.times.npy')) && ...
                ~isfile(fullfile(f, 'circles.times.npy'))
            continue
        end

        % load data
        spikeData = io.getEphysData(f);
        chanCoord = readNPY(fullfile(f, 'channels.localCoordinates.npy'));
        SC_depth = readNPY(fullfile(f, '_ss_recordings.scChannels.npy'));
        SC_top = chanCoord(SC_depth(1), 2);
        SC_SO = chanCoord(SC_depth(2), 2);
        results = cell(2,1);
        if isfile(fullfile(f, '_ss_sparseNoise.times.npy'))
            results{1} = io.getNoiseRFFits(f);
        end
        if isfile(fullfile(f, 'circles.times.npy'))
            results{2} = io.getCircleRFFits(f);
        end

        % select units with good RFs from noise or circles
        stimTypes = rf.selectRFStim(results, minEV, minPeak);
        
        % determine depth of each unit
        depths = spikeData.clusterDepths;

        % consider only units within SC and with good RF
        validUnits = find(depths(:,2) > 0 & ~isnan(stimTypes));
        
        % determine stimulus edges
        edges = NaN(2, 4);
        if ~isempty(results{1})
            edges(1,:) = results{1}.edges;
        end
        if ~isempty(results{2})
            gridX = results{2}.x;
            gridY = results{2}.y;
            gridW = median(diff(gridX));
            gridH = median(-diff(gridY));
            % edges: [left right top bottom] (above horizon: >0)
            edges(2,:) = [gridX(1)-0.5*gridW gridX(end)+0.5*gridW ...
                gridY(1)+0.5*gridH gridY(end)-0.5*gridH];
        end
        edges(1, [1 4]) = min(edges(:,[1 4]), [], 1);
        edges(1, [2 3]) = max(edges(:,[2 3]), [], 1);
        edges = edges(1,:);

        % determine colour for each unit
        colors = jet(500);
        [~, ~, colInds] = histcounts(depths(validUnits,1), ...
            linspace(0, SC_extent, 500));
        lines = {'-', ':'};

        figure
        tiledlayout(1, 1)
        nexttile
        hold on
        for k = 1:length(validUnits)
            unit = validUnits(k);
            res = results{stimTypes(unit)};
            if res.EV(unit) < minEV || res.peaks(unit) < minPeak
                continue
            end
            pars = res.gaussPars(unit,:);
            % ellipse at 2 STD (x and y), not rotated, not shifted
            x = pars(3) * cos(ellipse_x);
            y = pars(5) * sin(ellipse_x);
            % rotate and shift ellipse
            x_rot = pars(2) + x .* cos(pars(6)) - y .* sin(pars(6));
            y_rot = pars(4) + x .* sin(pars(6)) + y .* cos(pars(6));
            plot(x_rot, y_rot, lines{stimTypes(unit)}, ...
                'Color', colors(colInds(k),:), 'LineWidth', 2)
        end
        axis image
        axis(edges([1 2 4 3]))
        xlabel('Azimuth (visual degrees)')
        ylabel('Elevation (visual degrees)')

        nexttile('East')
        hold on
        imagesc(0, [0 SC_extent], (1:500)')
        plot([-0.5 0.5], [1 1] .* (SC_top - SC_SO), 'k', 'LineWidth', 2)
        colormap(colors)
        axis([-0.5 0.5 0 SC_extent])
        set(gca, "YAxisLocation", "right", "XTick", [], "Box", "off", ...
            "YDir", "reverse")
        ylabel('Depth from SC surface (in um)')

        sgtitle(sprintf('%s %s', name, date))
        
        saveas(gcf, fullfile(fPlots, sprintf('%s_%s_2D.jpg', name, date)));
        close gcf

        figure
        tiledlayout(1, 2)
        for dim = 1:2
            nexttile
            hold on
            for k = 1:length(validUnits)
                unit = validUnits(k);
                res = results{stimTypes(unit)};
                if res.EV(unit) < minEV || res.peaks(unit) < minPeak
                    continue
                end
                pars = res.gaussPars(unit,:);
                plot(pars(2*dim) + [-1 1] .* pars(2*dim + 1), ...
                    [1 1] .* depths(unit,1), ...
                    ['k' lines{stimTypes(unit)}], 'LineWidth', 2)
            end
            ylim([0 SC_extent])
            set(gca, 'YDir', 'reverse')
            if dim == 1
                plot(edges([1 2]), [1 1] .* SC_top - SC_SO, 'r')
                xlim(edges([1 2]))
                xlabel('Azimuth (visual degrees)')
                ylabel('Depth from SC surface (in um)')
            else
                plot(edges([4 3]), [1 1] .* SC_top - SC_SO, 'r')
                xlim(edges([4 3]))
                xlabel('Elevation (visual degrees)')
            end
        end

        sgtitle(sprintf('%s %s', name, date))
        
        saveas(gcf, fullfile(fPlots, sprintf('%s_%s_1D.jpg', name, date)));
        close gcf
    end
end