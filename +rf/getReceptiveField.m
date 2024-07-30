function receptiveFields = ...
    getReceptiveField(traces, traceTimes, ...
    stimFrames, stimTimes, RFtimesInFrames, ...
    lambda)

%GETRECEPTIVEFIELD Returns spatiotemporal receptive field.
%   [receptiveFields, explainedVariance, predictions, time] = ...
%    GETRECEPTIVEFIELD(traces, traceTimes, ...
%    stimFrames, stimTimes, RFtimesInFrames, ...
%    lambdas, crossFolds) calculates the linear RF of the neuron.
%
%   receptiveFields     [rows x cols x RFframes x RFtype x neuron]
%                       containing linear regression solution for x in Ax=B
%                       where A is stimulus [rows x cols x time] and B is
%                       calcium response, for each neuron and stimulus 
%                       model; ridge regression is performed on all data
%                       using the optimal lambda value found with
%                       cross-validation
%   explainedVariance   [neuron x lambdaStim x crossFold], each entry:
%                       explained variance for fitted RF for
%                       each neuron, lambda, and cross val. fold
%   predictions         [t x neuron], each column contains
%                       prediction based on RF for test 
%                       responses of specific neuron (using optimal
%                       lambda)
%   time                [t x 1]; time points for predictions
%
%   traces              [trTime x neuron]; calcium traces of neurons
%   traceTimes          [trTime x 1]; sample times of calcium traces
%   stimFrames          [time x rows x cols]; noise stimulus
%   stimTimes           [time x 1]; times of stimulus frames
%   RFtimesInFrames     [1 x RFframes]; frames of receptive field relative
%                       to stimulus frames
%   lambdas             [1 x lambda]; values of lambda
%   crossFolds          ind; number of cross val. folds

% generate toplitz matrix for stimulus
[stim, time, stimFrames, stimBin] = ...
    rf.makeStimToeplitz(stimFrames, stimTimes, RFtimesInFrames);

% get neural response
traceBin = median(diff(traceTimes));
numBins = round(stimBin / traceBin);
traces = smoothdata(traces, 1, 'movmean', numBins, 'omitnan');
traces = interp1(traceTimes, traces, time);
% z-score neural response
zTraces = (traces - nanmean(traces,1)) ./ nanstd(traces,0,1);

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
stim2 = (stim2 - nanmean(stim2(:))) ./ nanstd(stim2(:)); % normalise each column of stimulus matrix
clear sdesl

if isempty(lambda)
    lamStim = 0;
    lamMatrix_stim = [];
else
    % scale lamdas according to number of samples and number of predictors
    lamStim = sqrt(lambda .* size(stim,1) .* size(stim,2));

    % construct spatial smoothing lambda matrix
    lamMatrix_stim = rf.makeLambdaMatrix([size(stimFrames,2), size(stimFrames,3), ...
        length(RFtimesInFrames)], [1 1 0]);
    lamMatrix_stim = blkdiag(lamMatrix_stim, lamMatrix_stim);
end

% determine RFs using all data and optimal lambdas
y_train = gpuArray(padarray(zTraces(:,valid), size(lamMatrix_stim,1), 'post'));
x = stim2;

lms = lamMatrix_stim .* lamStim;

A = gpuArray([x; lms]);

receptiveFields = gather(A \ y_train);

receptiveFields = reshape(receptiveFields, size(stimFrames,2), ...
    size(stimFrames,3), length(RFtimesInFrames), 2, size(traces,2));