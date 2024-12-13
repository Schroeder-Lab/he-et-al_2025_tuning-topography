function [prediction, explainedVar] = predictFromRF(caTraces, t_ca, ...
    toeplitz, t_toeplitz, spatTempRF)
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

% get neural response sampled at stimulus frame times
tBin_ca = median(diff(t_ca));
tBin_stim = median(diff(t_toeplitz));
numBins = round(tBin_stim / tBin_ca);
caTraces = smoothdata(caTraces, 1, 'movmean', numBins, 'omitnan');
caTraces = interp1(t_ca, caTraces, t_toeplitz);
% z-score neural response
zTrace = (caTraces - mean(caTraces,1,'omitnan')) ./ std(caTraces,0,1,'omitnan');

% delete stim frames for which neuron has NaN
ind = isnan(zTrace);
toeplitz(ind,:) = [];
zTrace(ind,:) = [];

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

prediction = NaN(size(ind));
% multiply stimulus with RF to predict response
prediction(~ind) = stim * spatTempRF(:);
% get explained variance of prediction
explainedVar = 1 - sum((zTrace - prediction(~ind)).^2, 1) ./ sum(zTrace.^2, 1);