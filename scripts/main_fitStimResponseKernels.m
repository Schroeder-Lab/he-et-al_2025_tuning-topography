function main_fitStimResponseKernels(folders)
% Given the pre-processed calcium traces, fit response kernels, response
% amplitudes (per trial), and kernel shift (per stimulus, in case of bar 
% stimuli) to calcium responses to drifting gratings, static gratings, and
% moving bars.

%% Parameters
kernelLength = 10; % length of kernel in sec

%% Fit kernels
sets = {'neurons', 'boutons'};
stimTypes = {'gratingsDrifting', 'gratingsStatic', 'bars'};
for s = 1:2 % neurons and boutons
    subjDirs = dir(fullfile(folders.data, sets{s}, 'SS*'));
    for subj = 1:length(subjDirs) % animals
        name = subjDirs(subj).name;
        fprintf('%s\n', name)
        dateDirs = dir(fullfile(folders.data, sets{s}, name, '2*'));
        for dt = 1:length(dateDirs) %dates
            date = dateDirs(dt).name;
            fprintf('  %s\n', date)
            f = fullfile(folders.data, sets{s}, name, date);
            for k = 1:length(stimTypes)
                type = stimTypes{k};
                % ignore session if stimulus was not presented
                if ~isfile(fullfile(f, sprintf('_ss_%s.intervals.npy', type)))
                    continue
                end
                fprintf('    %s:', type)
                fPlots = fullfile(folders.plots, 'FittedKernels', type, ...
                    sets{s}, name, date);
                if ~isfolder(fPlots)
                    mkdir(fPlots)
                end
                
                % load data
                data = io.getCalciumData(f);
                time_traces = data.time;
                tr = data.traces;
                planes = data.planes;
                cellIDs = data.ids;
                delays = data.delays;
                
                data = io.getGratingInfo(f, type);
                time_stim = data.times;
                stimIDs = data.ids;
                stimRecording = data.interval;

                doShift = false;
                switch type
                    case 'gratingsDrifting'
                        stimDirs = data.directions;
                    case 'gratingsStatic'
                        stimDirs = data.orientations;
                    case 'bars'
                        stimDirs = data.directions;
                        doShift = true;
                end

                % ignore calcium traces before and after stimulus paradigm
                ind = time_traces >= time_stim(1,1)-1 & ...
                    time_traces <= time_stim(end,end)+1;
                t = time_traces(ind);
                tr = tr(ind,:);
                % subtract 8th percentile of each trace
                tr = tr - prctile(tr, 8, 1);

                % ignore "blank"/gray stimuli
                sampleRate = 1 / median(diff(time_traces));
                blankStim = find(isnan(stimDirs));
                indNoBlank = stimIDs ~= blankStim;
                stimsNoBlank = stimIDs(indNoBlank);
                % get stimulus onset times and duration in number of
                % samples
                stimOnsets = time_stim(indNoBlank,1);
                stimFrames = ceil(median(diff(time_stim,1,2)) * sampleRate);

                % length of kernel in number of samples
                kernelLengthSamples = round(kernelLength * sampleRate);
                % initialize results structure
                fitResults = struct('kernel', ...
                    repmat({NaN(kernelLengthSamples,1)}, size(tr,2), 1), ...
                    'amplitudes', NaN(size(time_stim,1)/length(stimDirs), ...
                    length(stimDirs)-length(blankStim)), ...
                    'lags', NaN(size(time_stim,1)/length(stimDirs), ...
                    length(stimDirs)-length(blankStim)), ...
                    'prediction', NaN(size(tr,1),1), ...
                    'pValue', NaN, 'R2', NaN);
                % loop over all ROIs
                for iCell = 1:size(tr,2)
                    if mod(iCell,20) == 0
                        fprintf(' %d', iCell)
                    end
                    % fit kernel, amplitudes, and kernel shifts
                    result = krnl.fitKernelIteratively(...
                        t + delays(planes(iCell)), tr(:,iCell), ...
                        stimOnsets, stimsNoBlank, kernelLengthSamples, ...
                        doShift, true, true, stimFrames);
                    fitResults(iCell) = result;
                    
                    % save and close plot
                    if result.pValue < 2
                        saveas(gcf, fullfile(fPlots, sprintf('Unit%03d.jpg', iCell)))
                    end
                    close(gcf)
                end
                fprintf('\n')
                % save results
                times.kernel = (0:kernelLengthSamples-1) ./ sampleRate;
                times.prediction = t;
                io.writeKernelFitResults(fitResults, times, f, type);
            end
        end
    end
end