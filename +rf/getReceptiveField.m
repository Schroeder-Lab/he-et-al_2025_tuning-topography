function receptiveFields = getReceptiveField(zTraces, toeplitz, ...
    stimSize, rfBins, lambda)
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

% clean up neural traces (delete times where all traces are NaN; if NaN 
% values < 10% in a neuron, exchange NaNs for 0; skip neurons that have only 
% NaN values)
validTimes = ~all(isnan(zTraces),2);
toeplitz(~validTimes,:) = [];
zTraces(~validTimes,:) = [];
% if NaN values < 10% in a neuron, exchange NaNs for 0
ind = any(isnan(zTraces),1) & sum(isnan(zTraces),1)/size(zTraces,1) <= 0.1;
if sum(ind) > 0
    zTraces(:,ind) = fillmissing(zTraces(:,ind),'constant',0);
end
% skip neurons that have only NaN values
validUnits = ~all(isnan(zTraces),1)';

% duplicate stimulus matrix to predict ON part (1st half) and OFF
% part (2nd half)
s = toeplitz;
s(toeplitz < 0) = 0;
stim = s;
s = toeplitz;
s(toeplitz > 0) = 0;
stim = [stim, s];
% normalise each column of stimulus matrix
stim = (stim - mean(stim(:),'omitnan')) ./ std(stim(:),'omitnan');

if isempty(lambda)
    lamStim = 0;
    lamMatrix_stim = [];
else
    % scale lamda according to number of samples and number of predictors
    lamStim = sqrt(lambda .* sum(validTimes) .* size(stim,2));

    % construct spatial smoothing lambda matrix
    lamMatrix_stim = rf.makeLambdaMatrix([stimSize, length(rfBins)], ...
        [1 1 0]);
    lamMatrix_stim = blkdiag(lamMatrix_stim, lamMatrix_stim);
end

% determine RFs
A = gpuArray([stim; lamMatrix_stim .* lamStim]);
y_train = gpuArray(padarray(zTraces(:,validUnits), ...
    size(lamMatrix_stim,1), 'post'));
receptiveFields = gather(A \ y_train);

receptiveFields = reshape(receptiveFields, ...
    [stimSize, length(rfBins), 2, size(zTraces,2)]);