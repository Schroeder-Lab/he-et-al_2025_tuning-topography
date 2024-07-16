%% Folders
getFolders;

%% Parameters
doShuffle = 1; % test whether calcium trace is better fit using a kernel
               % decribing stimulus response or just a baseline
doPlot = 1; % if 1, plot results for each cell (kernel fit and relationship
            % between nonvisual signal, stimulus response and baseline)

kernelLength = 10; %14; % length of kernel in sec (only for 'iterative' model)
basisLength = 0; %25; % length of cosine basis functions to model baseline in 
                  % sec (only for 'iterative' model)

% model evaluation
minR2 = 0.3;
minEV = 0;

%% Add paths
addpath(genpath(fullfile(folders.tools, 'npy-matlab')))
addpath(fullfile(folders.repo))

%% Fit kernels
sets = {'neurons', 'boutons'};
stimTypes = {'gratingsDrifting', 'gratingsStatic', 'bars'};
for s = 1:2
    subjDirs = dir(fullfile(folders.data, sets{s}, 'SS*'));
    for subj = 1:length(subjDirs)
        name = subjDirs(subj).name;
        dateDirs = dir(fullfile(folders.data, sets{s}, name, '2*'));
        for dt = 1:length(dateDirs)
            date = dateDirs(dt).name;
            f = fullfile(folders.data, sets{s}, name, date);
            for k = 1:length(stimTypes)
                type = stimTypes{k};
                if ~isfile(fullfile(f, sprintf('_ss_%s.intervals.npy', type)))
                    continue
                end
                fPlots = fullfile(folders.plots, type, sets{s}, name, date);
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
                stimDirs = data.directions;
                stimRecording = data.interval;

                % only consider data during stimuli
                ind = time_traces >= time_stim(1,1)-1 & ...
                    time_traces <= time_stim(end,end)+1;
                t = time_traces(ind);
                tr = tr(ind,:);
%                 % [upsample and time-shift traces]
%                 [t, tr] = traces.alignSampling(t, tr, delays, planes);
                % subtract 8th percentile
                tr = tr - prctile(tr, 8, 1);

                sampleRate = 1 / median(diff(time_traces));

                blankStim = find(isnan(stimDirs));
                indNoBlank = stimIDs ~= blankStim;
                stimsNoBlank = stimIDs(indNoBlank);
                stimOnsets = time_stim(indNoBlank,1);
                stimFrames = ceil(median(diff(time_stim,1,2)) * sampleRate);

                kernelLengthSamples = round(kernelLength * sampleRate);
                fitResults = struct('kernel', ...
                    repmat({NaN(kernelLengthSamples,1)}, size(tr,2), 1), ...
                    'amplitudes', NaN(size(time_stim,1)/length(stimDirs), ...
                    length(stimDirs)-length(blankStim)), ...
                    'prediction', NaN(size(tr,1),1), ...
                    'pValue', NaN, 'R2', NaN);
                for iCell = 1:size(tr,2)
                    result = krnl.fitKernelIteratively(...
                        t + delays(planes(iCell)), tr(:,iCell), ...
                        stimOnsets, stimsNoBlank, kernelLengthSamples, ...
                        true, true, stimFrames);
                    fitResults(iCell) = result;
                    
                    % save and close plot
                    if result.pValue < 2
                        saveas(gcf, fullfile(fPlots, sprintf('Unit%03d.jpg', iCell)))
                    end
                    close(gcf)
                end
                % save results
                times.kernel = (0:kernelLengthSamples-1) ./ sampleRate;
                times.prediction = t;
                io.writeKernelFitResults(fitResults, times, f, type);
            end
        end
    end
end