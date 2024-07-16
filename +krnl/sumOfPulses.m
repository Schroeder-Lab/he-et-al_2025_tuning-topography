function [amplitudes, kernel, prediction, R2] = ...
    sumOfPulses(trace, onsets, len, lenEssential)
% function [amplitudes, kernel] = SumOfPulses(trace, PulseOnsets, PulseLen, Pulse, RealAmps)
% splits the tracenal trace into a sum of pulses of a single estimated shape kernel
% of user-specified length PulseLen, occuring at user-specified times PulseOnsets, 
% with estimated amplitudes amplitudes.
% StimID is just for plotting purposes, and the last two arguments are optional, only for if you test it on synthetic data

if nargin < 4
    lenEssential = len;
end

trace = trace(:);
indNaN = find(isnan(trace));
traceLen = length(trace);

onsets = onsets(:);
nPulses = length(onsets);

BreakTol = 0.01;
MaxIter = 500;

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

% for the traceLen by nPulses Toeplitz Matrix, where trace ~ Toep*PulseShape
xt = repmat(1:len,mPulses,1);
yt = repmat(onsets(indValid),1,len) + repmat(0:(len-1),mPulses,1);
tUseMe = find(yt<=traceLen & ~ismember(yt, indNaN)); % this makes sure we don't write outside the matrix

% for the Shifted Pulse matrix, where trace ~ ShiftedPulses*PulseAmps
xs = repmat(1:mPulses,len,1);
ys = repmat((0:len-1)',1,mPulses) + repmat(onsets(indValid)',len,1);
sUseMe = find(ys<=traceLen & ~ismember(ys, indNaN)); % this makes sure we don't write outside the matrix

traceNoNaN = trace;
traceNoNaN(indNaN) = 0;
for i=1:MaxIter
    % make the toeplitz matrix based on current estimated amplitudes
    st = repmat(amplitudes,1,len);
    Toep = sparse(yt(tUseMe), xt(tUseMe), st(tUseMe), traceLen, len);

    % simple linear regression
    b = Toep\traceNoNaN;
    kernel = b(1:len);
    
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
    b = ShiftedPulses\traceNoNaN;
    amplitudes = b(1:mPulses);
    
    % break if tolerance achieved
    if norm(amplitudes-Oldamplitudes)<BreakTol; break; end
end
amps = NaN(nPulses,1);
amps(indValid) = amplitudes;
amplitudes = amps;
prediction = ShiftedPulses * b;

R2 = 1 - sum((traceNoNaN - prediction).^2) / ...
    sum((traceNoNaN - mean(traceNoNaN)).^2);
% adjR2 = 1 - (sum((traceNoNaN - prediction).^2) / ...
%     (traceLen - length(indNaN) - len - mPulses)) / ...
%     (sum((traceNoNaN - mean(traceNoNaN)).^2) / (traceLen - length(indNaN) - 1));