function fitResults = ...
    fitKernelIteratively(time, trace, stimOnsets, stimIDs, kernelLength, ...
    doPlot, doShuffle, stimLength)

% time                  [time x 1], time of trace samples
% trace                 [time x 1], trace of one neuron
% stimOnsets            [t x 1], times of stimulus onsets, blanks not
%                       included
% stimIDs               [t x 1], stimulus IDs
% kernelLength          int; number of samples defining duration of kernel
% doPlot                if 1, data and kernels are plotted
% doShuffle             if 1, calcium trace is shifted randomly to estimate
%                       null hypothesis of kernels
% stimLength            int, stimulus duration in samples (time *
%                       samplingRate); only important if kernelLength much
%                       longer than stimLength

testRep = 500;
sigThreshold = 0.05;
batch = 50;
numBatches = ceil(testRep/batch);

if nargin < 5
    doPlot = 0;
end
if nargin < 6
    doShuffle = 0;
end
if nargin < 7
    stimLength = kernelLength;
end

% get stimulus onsets in frames
stimOnsetFrames = find(histcounts(stimOnsets, time));

[alphas, kernel, prediction, R2] = ...
    krnl.sumOfPulses(trace, stimOnsetFrames, kernelLength, stimLength);

stimUnique = unique(stimIDs);
numTrials = ceil(length(stimOnsets) / length(stimUnique));
if ~isempty(alphas)
    alphaMat = NaN(numTrials, length(stimUnique));
    % [trials x stimuli]
    for iStim = 1:length(stimUnique)
        alphaMat(:,iStim) = alphas(stimIDs == stimUnique(iStim));
    end
else
    alphaMat = [];
end

fitResults.kernel = kernel(:);
fitResults.amplitudes = alphaMat;
fitResults.prediction = prediction;
fitResults.pValue = NaN;
fitResults.R2 = R2;

% shift calcium trace and redo kernel to test whether cell has stimulus
% response
if doShuffle == 1
    shuffled = mod(bsxfun(@plus, randi(length(trace), 1, testRep), ...
        (0:length(trace)-1)'), length(trace));
    shuffled(shuffled==0) = length(trace);
    shuffled = trace(shuffled);
    R2Shuffled = NaN(testRep,1);
    for iBatch = 1:numBatches
        n = min(batch, testRep-(iBatch-1)*batch);
        r2Sh = zeros(1, n);
        sh = shuffled(:,(1:n)+(iBatch-1)*batch);
        parfor iRep = 1:n
            [~,~,~,r] = krnl.sumOfPulses(sh(:,iRep), stimOnsetFrames, ...
                kernelLength, stimLength);
            r2Sh(iRep) = r;
        end
        R2Shuffled((1:n)+(iBatch-1)*batch) = r2Sh;
        if iBatch < numBatches && sum(R2Shuffled > R2) / testRep > sigThreshold
            fitResults.pValue = 2;
            break
        end
    end

    if ~any(isnan(R2Shuffled))
        fitResults.pValue = sum(R2Shuffled > R2) / testRep;
    end
end

if doPlot == 1
    binSize = median(diff(time));

    figure('Position', [1921 1 1920 1123]);
    % plots some stuff....
    subplot(3,2,1)
    plot((0:kernelLength-1) .* binSize, kernel)
    hold on
    plot((stimLength-1) .* binSize, kernel(stimLength), 'ro')
    xlim([0 kernelLength-1] .* binSize)
    title('Response waveform');
    xlabel('time (s)')
    subplot(3,1,2), cla
    stem(stimOnsets, alphas, 'r', 'filled')
    label = {'stim resp'};
    legend(label);
    title('Estimated amplitudes')
    ax1=gca;
    subplot(3,1,3); cla
    plot(time, trace, 'k'); hold on
    plot(time, prediction, 'r');
    label = {'original', 'reconstructed'};
    plot([stimOnsets stimOnsets]', repmat(ylim,length(stimOnsets),1)', 'c:');
    legend(label);
    xlabel('Samples')
    title(sprintf('Original and fitted traces (R2: %.2f, p = %.3f)', ...
        R2, fitResults.pValue))
    ax2=gca;
    linkaxes([ax1 ax2],'x');
    xlim(time([1 end]))
end