function [prediction, explainedVar] = predictFromRF(trace, traceTimes, ...
    stimFrames, stimTimes, RFtimesInFrames, spatTempRF)
%PREDICTFROMRF   Predict response trace from stimulus and receptive field.

% INPUTS
% trace             [t], traces of unit
% traceTime         [t], sample times of traces
% stimFrames        [t_st x rows x cols]; noise stimulus
% stimTimes         [t_st]; times of stimulus frames
% RFtimesInFrames   [1 x RFframes]; frames of receptive field relative
% spatTempRF        [rows x columns x t_rf x ON/OFF], 2D masks of fitted ON
%                   and OFF subfields

% OUTPUTS
% prediction        [t], predicted trace
% explainedVar      double, explained variance

% generate toplitz matrix for stimulus
[stim, time, ~, stimBin] = ...
    rf.makeStimToeplitz(stimFrames, stimTimes, RFtimesInFrames);

% get neural response sampled at stimulus frame times
traceBin = median(diff(traceTimes));
numBins = round(stimBin / traceBin);
trace = smoothdata(trace, 1, 'movmean', numBins, 'omitnan');
trace = interp1(traceTimes, trace, time);
% z-score neural response
zTrace = (trace - mean(trace,1,'omitnan')) ./ std(trace,0,1,'omitnan');

% delete stim frames for which neuron has NaN
ind = isnan(zTrace);
stim(ind,:) = [];
zTrace(ind,:) = [];

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

prediction = NaN(size(ind));
% multiply stimulus with RF to predict response
prediction(~ind) = stim2 * spatTempRF(:);
% get explained variance of prediction
explainedVar = 1 - sum((zTrace - prediction(~ind)).^2, 1) ./ sum(zTrace.^2, 1);