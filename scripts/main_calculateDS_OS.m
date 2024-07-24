% Given the fitted response kernels (from main_fitStimResponseKernels.m),
% use vector averaging to determine preferred direction and orientation 
% (direction of final vector) as well as direction and orientation 
% selectivity (length of final vector). Use permutation test to determine
% significance of direction and orientation selectivity.

%% Folders
getFolders;

%% Parameters
numShuffles = 1000; % number of permutations for testing significance of 
                    % direction and orientation selectivity
maxP = 0.05; % p-value threshold for response kernel and 
             % direction/orientation selectivity
minR2 = 0; % threshold for explained variance for response kernel

% colors for plotting
cols = [1 0 0; 1 0.8 0.8];
phaseCols = lines(6);

%% Add paths
addpath(genpath(fullfile(folders.tools, 'npy-matlab')))
addpath(fullfile(folders.repo))

%% Plot kernel, direction and orientation tuning curves
sets = {'neurons', 'boutons'};
stimTypes = {'gratingsDrifting', 'bars', 'gratingsStatic'};
for s = 1:2 % neurons and boutons
    subjDirs = dir(fullfile(folders.data, sets{s}, 'SS*'));
    for subj = 1:length(subjDirs) % animals
        name = subjDirs(subj).name;
        dateDirs = dir(fullfile(folders.data, sets{s}, name, '2*'));
        for dt = 1:length(dateDirs) %dates
            date = dateDirs(dt).name;
            f = fullfile(folders.data, sets{s}, name, date);
            for k = 1:length(stimTypes)
                type = stimTypes{k};
                % ignore session if stimulus was not presented
                if ~isfile(fullfile(f, sprintf('_ss_%s.intervals.npy', type)))
                    continue
                end
                fPlots = fullfile(folders.plots, 'TuningCurves', type, ...
                    sets{s}, name, date);
                if ~isfolder(fPlots)
                    mkdir(fPlots)
                end
                
                % load data                
                data = io.getStimResponseFits(f, type);
                kernel = data.kernel;
                time_kernel = data.time_kernel;
                amplitudes = data.amplitudes;
                p_kernel = data.pValue;
                r2_kernel = data.R2;

                data = io.getGratingInfo(f, type);
                switch type
                    case {'gratingsDrifting', 'bars'}
                        stimDirs = data.directions;
                        % ignore blank stimulus
                        stimDirs(isnan(stimDirs)) = [];
                        % determine values needed for plotting
                        stimDist = median(diff(stimDirs));
                        stimDirsCirc = [stimDirs; 360];
                        stimOrisCirc = stimDirs(1:length(stimDirs)/2+1);
                    case 'gratingsStatic'
                        stimDirs = data.orientations;
                        stimPhases = data.phases;
                        % ignore blank stimulus
                        indNaN = isnan(stimDirs);
                        stimDirs(indNaN) = [];
                        stimPhases(indNaN) = [];
                        % determine values needed for plotting
                        phases = unique(stimPhases);
                        ind0 = find(stimDirs == 0);
                        stimOrisCirc = [stimDirs; ...
                            ones(length(ind0),1) .* 180];
                        orientations = unique(stimOrisCirc);
                        stimDist = median(diff(orientations));
                        stimPhasesCirc = [stimPhases; stimPhases(ind0)];
                end
                stimTimes = data.times;
                stimDur = median(diff(stimTimes,1,2));

                % initialize results structures
                dirTuning = struct('preference', ...
                    repmat({NaN(1,1)}, size(kernel,2), 1), ...
                    'selectivity', NaN, 'pValue', NaN, 'responseSign', NaN);
                oriTuning = struct('preference', ...
                    repmat({NaN(1,1)}, size(kernel,2), 1), ...
                    'selectivity', NaN, 'pValue', NaN, 'responseSign', NaN);
                indResponsive = find(p_kernel < maxP & r2_kernel > minR2)';
                % loop over all ROIs
                for iCell = indResponsive
                    % determine direction/orientation preference,
                    % selectivity, and signficance (response sign: -1 if
                    % most responses are negative)
                    [drct, ori] = tuning.determineTuning(amplitudes(:,:,iCell), ...
                        stimDirs, numShuffles);
                    dirTuning(iCell) = drct;
                    oriTuning(iCell) = ori;
                    maxi = max(amplitudes(:,:,iCell),[],'all');
                    mini = min(amplitudes(:,:,iCell),[],'all');

                    % plot results
                    % replicate responses at 0 degrees to plot at 360 or
                    % 180 degrees
                    switch type
                        case {'gratingsDrifting', 'bars'}
                            ampsDir = amplitudes(:, [1:end,1], iCell) .* drct.responseSign;
                            ampsOri = [amplitudes(:,1:length(stimDirs)/2, iCell); ...
                                amplitudes(:,length(stimDirs)/2+1:end, iCell)] .* ...
                                ori.responseSign;
                            ampsOri = ampsOri(:, [1:end 1]);
                            figure('Position', [100 380 1740 420]);
                            tiledlayout(1,3)
                        case 'gratingsStatic'
                            ampsOri = amplitudes(:, [1:end,find(stimDirs==0)'], ...
                                iCell) .* ori.responseSign;
                            figure('Position', [100 380 1200 420]);
                            tiledlayout(1,2)
                    end
                    nexttile % response kernel
                    plot(time_kernel, kernel(:,iCell), 'k', 'LineWidth', 2)
                    hold on
                    plot([1 1].*stimDur, [min([0 min(kernel(:,iCell))]) 1], ...
                        'k:')
                    xlim(time_kernel([1 end]))
                    title(sprintf('Kernel (R2 = %.2f, p = %.3f)', ...
                        r2_kernel(iCell), p_kernel(iCell)))
                    xlabel('Time (s)')
                    set(gca, 'box', 'off')
                    
                    % direction tuning curve (not if grating was static)
                    if any(strcmp(type, stimTypes(1:2)))
                        nexttile
                        c = cols(1,:);
                        if drct.pValue > maxP
                            c = cols(2,:);
                        end
                        plot(stimDirsCirc' + randn(size(ampsDir)) .* (stimDist/10), ...
                            ampsDir, '.', 'MarkerSize', 10, ...
                            'Color', [1 1 1].*0.5);
                        hold on
                        plot(stimDirsCirc, mean(ampsDir,1, 'omitnan'), ...
                            'k.-', 'MarkerSize', 30, 'LineWidth', 1);
                        plot(drct.preference, maxi, 'v', 'MarkerSize', 8, ...
                            'MarkerFaceColor', c, 'MarkerEdgeColor', 'none');
                        xlim([-10 370])
                        ylim([min([0 1.1*mini]) 1.1*maxi])
                        xlabel('Direction (deg)')
                        title(sprintf('DS = %.3f, p = %.3f', drct.selectivity, ...
                            drct.pValue))
                        set(gca, 'box', 'off', 'XTick', stimDirsCirc(1:3:end))
                    end

                    nexttile % orientation tuning curve
                    c = cols(1,:);
                    h = [0 0 0];
                    if ori.pValue > maxP
                        c = cols(2,:);
                    end
                    if strcmp(type, stimTypes{3}) % for static gratings
                        hold on
                        % color code responses for difference phases
                        for ph = 1:length(phases)
                            indSt = stimPhasesCirc == phases(ph);
                            hs = plot(stimOrisCirc(indSt)' + ...
                                randn(size(ampsOri(:,indSt))) .* (stimDist/10), ...
                                ampsOri(:,indSt), '.', 'MarkerSize', 10, ...
                                'Color', phaseCols(ph,:));
                        end
                        h(1) = hs(1);
                        % average responses for same orientation across
                        % phases
                        meanAmps = NaN(length(orientations),1);
                        for st = 1:length(orientations)
                            indSt = stimOrisCirc == orientations(st);
                            meanAmps(st) = mean(ampsOri(:,indSt), 'all', 'omitnan');
                        end
                        h(2) = plot(orientations, meanAmps, ...
                            'k.-', 'MarkerSize', 30, 'LineWidth', 1);
                        h(3) = plot(ori.preference, maxi, 'v', 'MarkerSize', 8, ...
                            'MarkerFaceColor', c, 'MarkerEdgeColor', 'none');
                        legend(h, 'single trials per phase','median','preference', ...
                            'Location', 'bestoutside')
                        set(gca, 'box', 'off', 'XTick', orientations)
                    else
                        hs = plot(stimOrisCirc' + randn(size(ampsOri)) .* (stimDist/10), ...
                            ampsOri, '.', 'MarkerSize', 10, ...
                            'Color', [1 1 1].*0.5);
                        h(1) = hs(1);
                        hold on
                        h(2) = plot(stimOrisCirc, mean(ampsOri,1, 'omitnan'), ...
                            'k.-', 'MarkerSize', 30, 'LineWidth', 1);
                        h(3) = plot(ori.preference, maxi, 'v', 'MarkerSize', 8, ...
                            'MarkerFaceColor', c, 'MarkerEdgeColor', 'none');
                        legend(h, 'single trials','median','preference', ...
                            'Location', 'bestoutside')
                        set(gca, 'box', 'off', 'XTick', stimOrisCirc(1:2:end))
                    end
                    xlim([-10 190])
                    ylim([min([0 1.1*mini]) 1.1*maxi])
                    xlabel('Orientation (deg)')
                    title(sprintf('OS = %.3f, p = %.3f', ori.selectivity, ...
                        ori.pValue))
                    
                    % save and close plot
                    saveas(gcf, fullfile(fPlots, sprintf('Unit%03d.jpg', iCell)))
                    close(gcf)
                end
                % save data
                io.writeTuningResults(dirTuning, oriTuning, f, type);
            end
        end
    end
end