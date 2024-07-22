function fitResults = ...
    fitKernelIteratively(time, trace, stimOnsets, stimIDs, kernelLength, ...
    doShift, doPlot, doShuffle, stimLength)
%FITKERNELITERATIVELY   Fit stimulus responses with kernel, amplitudes and
%shifts, then perform shift test, and plot results.

% INPUTS
% time                  [time x 1], time of trace samples
% trace                 [time x 1], trace of one neuron
% stimOnsets            [t x 1], times of stimulus onsets, blanks not
%                       included
% stimIDs               [t x 1], stimulus IDs
% kernelLength          int; number of samples defining duration of kernel
% doShift               calculate optimal shift of kernel relative to
%                       stimulus onset; used for moving bars that enter RF
%                       at different latencies
% doPlot                if 1, data and kernels are plotted
% doShuffle             if 1, calcium trace is shifted randomly to estimate
%                       null distribution of explained variance
% stimLength            int, stimulus duration in samples (time *
%                       samplingRate); only important if kernelLength much
%                       longer than stimLength

% OUTPUTS
% fitResults
%   .kernel             [k x 1], shape of kernel, normalized to max = 1
%   .amplitudes         [rep x stim], kernel amplitudes for each repetition
%                       (rep) of each stimulus
%   .lags               [rep x stim], kernel lag for each repetition of
%                       each stimulus
%   .prediction         [time x 1], predicted calcium trace
%   .pValue             [1], p-value of fit based on shift-test; set to 2
%                       if more than 5% of current set of shifted data
%                       resulted in better fit than non-shifted data
%   .R2                 [1], explained variance of fit

%% Parameters
testRep = 500; % generate 500 trace shifts for significance test
sigThreshold = 0.05; % significance threshold to stop fit for further shifted traces
batch = 50; % number of shifts to test in each round
numBatches = ceil(testRep/batch);

%% Defaults
if nargin < 5
    doShift = false;
end
if nargin < 6
    doPlot = 0;
end
if nargin < 7
    doShuffle = 0;
end
if nargin < 8
    stimLength = kernelLength;
end

%% Run the fit
% get stimulus onsets in samples
stimOnsetFrames = find(histcounts(stimOnsets, time));
% fit data
[alphas, kernel, lags, prediction, R2] = ...
    krnl.sumOfPulses(trace, stimOnsetFrames, doShift, stimIDs, ...
    kernelLength, stimLength);
% reshape vector of amplitudes and lags
stimUnique = unique(stimIDs);
numTrials = ceil(length(stimOnsets) / length(stimUnique));
if ~isempty(alphas)
    alphaMat = NaN(numTrials, length(stimUnique));
    % [trials x stimuli]
    lagMat = NaN(numTrials, length(stimUnique));
    for iStim = 1:length(stimUnique)
        alphaMat(:,iStim) = alphas(stimIDs == stimUnique(iStim));
        lagMat(:,iStim) = lags(stimIDs == stimUnique(iStim));
    end
else
    alphaMat = [];
    lagMat = [];
end

fitResults.kernel = kernel(:);
fitResults.amplitudes = alphaMat;
fitResults.lags = lagMat;
fitResults.prediction = prediction;
fitResults.pValue = NaN;
fitResults.R2 = R2;

%% Run the shift-test
% shift calcium trace and redo kernel to test whether cell has stimulus
% response
if doShuffle == 1
    % create random shifts
    shuffled = mod(bsxfun(@plus, randi(length(trace), 1, testRep), ...
        (0:length(trace)-1)'), length(trace));
    shuffled(shuffled==0) = length(trace);
    shuffled = trace(shuffled);
    % fit shifted data and determine explained variance for each shift
    R2Shuffled = NaN(testRep,1);
    for iBatch = 1:numBatches
        n = min(batch, testRep-(iBatch-1)*batch);
        r2Sh = zeros(1, n);
        sh = shuffled(:,(1:n)+(iBatch-1)*batch);
        parfor iRep = 1:n
            [~,~,~,~,r] = krnl.sumOfPulses(sh(:,iRep), stimOnsetFrames, ...
                doShift, stimIDs, kernelLength, stimLength);
            r2Sh(iRep) = r;
        end
        R2Shuffled((1:n)+(iBatch-1)*batch) = r2Sh;
        % stop loop if max. p-value is reached
        if iBatch < numBatches && sum(R2Shuffled > R2) / testRep > sigThreshold
            fitResults.pValue = 2;
            break
        end
    end
    % update p-value
    if ~any(isnan(R2Shuffled))
        fitResults.pValue = sum(R2Shuffled > R2) / testRep;
    end
end

%% Plot the results
if doPlot == 1
    binSize = median(diff(time));

    figure('Position', [1921 1 1920 1123]);
    % plots kernel
    subplot(3,2,1)
    plot((0:kernelLength-1) .* binSize, kernel)
    hold on
    plot([1 1] .* (stimLength-1) .* binSize, [0 1], 'r')
    legend('kernel','stim offset')
    xlim([0 kernelLength-1] .* binSize)
    title('Response waveform');
    xlabel('time (s)')
    % plot amplitudes
    subplot(3,1,2), cla
    stem(stimOnsets, alphas, 'r', 'filled')
    label = {'stim resp'};
    legend(label);
    title('Estimated amplitudes')
    ax1=gca;
    % plot predicted calcium trace
    subplot(3,1,3); cla
    plot(time, trace, 'k'); hold on
    plot(time, prediction, 'r');
    label = {'original', 'reconstructed'};
    plot([stimOnsets stimOnsets]', repmat(ylim,length(stimOnsets),1)', 'c:');
    legend(label);
    xlabel('Time (s)')
    title(sprintf('Original and fitted traces (R2: %.2f, p = %.3f)', ...
        R2, fitResults.pValue))
    ax2=gca;
    linkaxes([ax1 ax2],'x');
    xlim(time([1 end]))
end