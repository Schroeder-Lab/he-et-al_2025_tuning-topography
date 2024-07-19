function [amplitudes, kernel, lags, prediction, R2] = ...
    sumOfPulses(trace, onsets, doShift, stimIDs, len, lenEssential)
% function [amplitudes, kernel] = SumOfPulses(trace, PulseOnsets, PulseLen, Pulse, RealAmps)
% splits the tracenal trace into a sum of pulses of a single estimated shape kernel
% of user-specified length PulseLen, occuring at user-specified times PulseOnsets, 
% with estimated amplitudes amplitudes.
% StimID is just for plotting purposes, and the last two arguments are optional, only for if you test it on synthetic data

if nargin < 6
    lenEssential = len;
end

maxShift = round(0.4 * lenEssential);
kernelThreshold = 0.1;
BreakTol = 0.01;
MaxIter = 500;

trace = trace(:);
indNaN = find(isnan(trace));
traceLen = length(trace);

onsets = onsets(:);
nPulses = length(onsets);

% The algorithm works by iterating two steps: computing the amplitudes
% then computing the pulse shape. Both of these are done by a backslash
% operation of a sparse matrix. To do this quickly, we first make up some 
% matrices that will help build it.

% for the Shifted Pulse matrix, where trace ~ ShiftedPulses * PulseAmps
xs = repmat(1:nPulses, lenEssential, 1); % column indices for pulse IDs
ys = repmat((0:lenEssential-1)', 1, nPulses) + repmat(onsets(:)', ...
    lenEssential, 1); % time indices for pulses

sUseMe = find(ys<=traceLen & ~ismember(ys, indNaN)); % this makes sure we don't write outside the matrix

ss = ones(nPulses, lenEssential);
ShiftedPulses = sparse(ys(sUseMe), xs(sUseMe), ss(sUseMe), ...
    traceLen, nPulses); % [time x pulse], 1 if pulse is occurring
indValid = sum(ShiftedPulses,1) >= .7*lenEssential;
mPulses = sum(indValid);

% initial guess for the amplitudes: all ones.
amplitudes = ones(mPulses,1);
% initial guess for lags: all zeros.
stims = unique(stimIDs);
lagsTrials = zeros(1,mPulses);

% for the traceLen by nPulses Toeplitz Matrix, where trace ~ Toep*PulseShape
xt = repmat(1:len,mPulses,1);
yt = repmat(onsets(indValid),1,len) + repmat(0:(len-1),mPulses,1);
tUseMe = find(yt<=traceLen & ~ismember(yt, indNaN)); % this makes sure we don't write outside the matrix

% for the Shifted Pulse matrix, where trace ~ ShiftedPulses*PulseAmps
xs = repmat(1:mPulses,len,1);
ys = repmat(onsets(indValid)',len,1) + repmat((0:len-1)',1,mPulses);
sUseMe = find(ys<=traceLen & ~ismember(ys, indNaN)); % this makes sure we don't write outside the matrix

traceNoNaN = trace;
traceNoNaN(indNaN) = 0;

ys_original = ys; 
yt_original = yt;
for i=1:MaxIter
    yt = yt_original + lagsTrials';
    tUseMe = find(yt<=traceLen & ~ismember(yt, indNaN) & yt>0); % this makes sure we don't write outside the matrix

    % make the toeplitz matrix based on current estimated amplitudes
    st = repmat(amplitudes,1,len);
    Toep = sparse(yt(tUseMe), xt(tUseMe), st(tUseMe), traceLen, len);

    % simple linear regression
    kernel = Toep \ traceNoNaN;
    
    % normalize to have peak one
    kernel = kernel / max(abs(kernel));
    % make positive going
    if max(kernel) < abs(min(kernel))
        kernel = -kernel;
    end

    % construct a matrix of shifted pulses
    ss = repmat(kernel,mPulses,1);
    ShiftedPulses = sparse(ys(sUseMe), xs(sUseMe), ss(sUseMe), ...
        traceLen, mPulses);
    
    % now estimate the amplitude of each pulse
    Oldamplitudes = amplitudes;
    
    % simple linear regression
    amplitudes = ShiftedPulses\traceNoNaN;
    
    % break if tolerance achieved
    if norm(amplitudes-Oldamplitudes)<BreakTol; break; end
    
    % CONTINUE HERE-------------------------------------------------
    % todo: is zscoring better? is dividing by max better? don't set lag to
    % 0 if it was larger than maxShift. instead set to median(lags)?
    % consider maxShift because lags are ADDED!
    % determine lags
    if doShift
        oldLags = lagsTrials;
        lagsTrials = NaN(1,mPulses);

        indUse = ys > 0 & ys <= traceLen;
        yUse = ys;
        yUse(~indUse) = 1;

        indMax = find(kernel == 1);
        indThr = find(kernel(1:indMax) < kernelThreshold, 1, 'last');
        if indThr > 2
            kernel = kernel([indThr-1:end, 1:indThr-2]);
        end
        for s = 1:length(stims)
            trials = stimIDs(indValid) == stims(s);
            pred = ShiftedPulses(:,~trials) * amplitudes(~trials);
            tr = trace - pred;
            pulses = tr(yUse);
            pulses(~indUse) = NaN;
            pulse = mean(pulses(:,trials), 2, 'omitnan');
            % pulse = pulse ./ max(abs(pulse));
            pulseModel = mean(kernel * amplitudes(trials)',2);
            % pulseModel = pulseModel ./ max(abs(pulseModel));
            [corrs, l] = xcorr(pulse, pulseModel, maxShift, 'unbiased');
            [~,indCorr] = max(corrs);
            lagsTrials(trials) = l(indCorr);

            % figure,plot((1:len)-l(indCorr),pulse),hold on,plot(pulseModel)
        end
        lagsTrials = oldLags + lagsTrials;
        lagsTrials = lagsTrials - max(lagsTrials);
        lagsTrials(lagsTrials < -2*maxShift) = round(median(lagsTrials(lagsTrials > -2*maxShift)));
        %update ys matrix based on lags
        ys = ys_original + lagsTrials;
        sUseMe = find(ys<=traceLen & ~ismember(ys, indNaN) & ys>0); %find(ys<=SigLen & ys>0)
    end
end
prediction = ShiftedPulses * amplitudes;
amps = NaN(nPulses,1);
amps(indValid) = amplitudes;
amplitudes = amps;
lags = zeros(nPulses,1);
lags(indValid) = lagsTrials';

R2 = 1 - sum((traceNoNaN - prediction).^2) / ...
    sum((traceNoNaN - mean(traceNoNaN)).^2);
% adjR2 = 1 - (sum((traceNoNaN - prediction).^2) / ...
%     (traceLen - length(indNaN) - len - mPulses)) / ...
%     (sum((traceNoNaN - mean(traceNoNaN)).^2) / (traceLen - length(indNaN) - 1));