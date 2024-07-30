function receptiveFields = getReceptiveField(traces, traceTimes, ...
    stimFrames, stimTimes, RFtimesInFrames, lambda)
%GETRECEPTIVEFIELD   Return response-triggered spatiotemporal receptive field.

% INPUTS
% traces              [t x ROIs]; calcium traces of units
% traceTimes          [t]; sample times of calcium traces
% stimFrames          [t_st x rows x cols]; noise stimulus
% stimTimes           [t_st]; times of stimulus frames
% RFtimesInFrames     [1 x RFframes]; frames of receptive field relative
%                     to response
% lambda              [1 x lambda]; regularisation lambda

% OUTPUTS
% receptiveFields     [rows x cols x RFframes x RFtype x ROIs]
%                     containing linear regression solution for x in Ax=B
%                     where A is stimulus [rows x cols x time] and B is
%                     calcium response, for each unit and ON/OFF response; 
%                     ridge regression is performed using lambda

% generate toplitz matrix for stimulus
[stim, time, stimFrames, stimBin] = ...
    rf.makeStimToeplitz(stimFrames, stimTimes, RFtimesInFrames);

% get neural response
traceBin = median(diff(traceTimes));
numBins = round(stimBin / traceBin);
traces = smoothdata(traces, 1, 'movmean', numBins, 'omitnan');
traces = interp1(traceTimes, traces, time);
% z-score neural response
zTraces = (traces - mean(traces,1,'omitnan')) ./ std(traces,0,1,'omitnan');

% delete stim frames for which all neurons have NaN
ind = all(isnan(zTraces),2);
stim(ind,:) = [];
zTraces(ind,:) = [];
% if NaN values < 5% in a neuron, exchange NaNs for 0
ind = any(isnan(zTraces),1) & sum(isnan(zTraces),1)/size(zTraces,1) <= 0.05;
if sum(ind) > 0
    zTraces(:,ind) = fillmissing(zTraces(:,ind),'constant',0);
end
% skip neurons that have only NaN values
valid = ~all(isnan(zTraces),1)';

% duplicate stimulus matrix to predict ON part (1st half) and OFF
% part (2nd half)
s = stim;
s(stim < 0) = 0;
stim2 = s;
s = stim;
s(stim > 0) = 0;
stim2 = [stim2, s];
% normalise each column of stimulus matrix
stim2 = (stim2 - mean(stim2(:),'omitnan')) ./ std(stim2(:),'omitnan');
clear sdesl

if isempty(lambda)
    lamStim = 0;
    lamMatrix_stim = [];
else
    % scale lamda according to number of samples and number of predictors
    lamStim = sqrt(lambda .* size(stim,1) .* size(stim,2));

    % construct spatial smoothing lambda matrix
    lamMatrix_stim = rf.makeLambdaMatrix([size(stimFrames,2), size(stimFrames,3), ...
        length(RFtimesInFrames)], [1 1 0]);
    lamMatrix_stim = blkdiag(lamMatrix_stim, lamMatrix_stim);
end

% determine RFs
y_train = gpuArray(padarray(zTraces(:,valid), size(lamMatrix_stim,1), 'post'));
x = stim2;

lms = lamMatrix_stim .* lamStim;

A = gpuArray([x; lms]);

receptiveFields = gather(A \ y_train);

receptiveFields = reshape(receptiveFields, size(stimFrames,2), ...
    size(stimFrames,3), length(RFtimesInFrames), 2, size(traces,2));