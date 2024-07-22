%% Folders
getFolders;

%% Parameters
numShuffles = 1000;

%% Add paths
addpath(genpath(fullfile(folders.tools, 'npy-matlab')))
addpath(fullfile(folders.repo))

%% Plot kernel, direction and orientation tuning curves
sets = {'neurons', 'boutons'};
stimTypes = {'gratingsDrifting', 'bars'};
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
                stimDirs = data.directions;

                % determine values needed for plotting
                stimDist = median(diff(stimDirs));

                % initialize results structures
                dirTuning = struct('preference', ...
                    repmat({NaN(1,1)}, size(kernel,2), 1), ...
                    'selectivity', NaN, 'pValue', NaN);
                oriTuning = struct('preference', ...
                    repmat({NaN(1,1)}, size(kernel,2), 1), ...
                    'selectivity', NaN, 'pValue', NaN);
                % loop over all ROIs
                for iCell = 1:size(tr,2)
                    [dir, ori] = tuning.determineTuning(amplitudes, ...
                        stimDirs, numShuffles);
                    dirTuning(iCell) = dir;
                    oriTuning(iCell) = ori;

                    figure('Position', [680 50 560 946]);
                    tiledlayout(3,1)
                    nexttile
                    plot(time_kernel, kernel, 'k')
                    xlim(time_kernel([1 end]))
                    title(sprtinf('Kernel (R2 = %.2f, p = %.3f)', ...
                        r2_kernel, p_kernel))
                    xlabel('Time (s)')
                    
                    % TODO: mark preferred direction
                    nexttile
                    plot(stimDirs + randn(size(amplitudes)) .* (stimDist/10), ...
                        amplitudes, '.', 'MarkerSize', 10, ...
                        'MarkerFaceColor', [1 1 1].*0.5, ...
                        'MarkerEdgeColor', 'none')
                    hold on
                    plot(stimDirs, median(amplitudes,1, 'omitnan'), ...
                        'k.--', 'MarkerSize', 15, 'LineWidth', 2)
                    xlabel('Direction (deg)')
                    title(sprintf('DS = %.3f, p = %.3f', dir.selectivity, ...
                        dir.pValue))
                    
                    % save and close plot
                    if result.pValue < 2
                        saveas(gcf, fullfile(fPlots, sprintf('Unit%03d.jpg', iCell)))
                    end
                    close(gcf)
                end
            end
        end
    end
end