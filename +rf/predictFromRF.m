function [prediction, explainedVar] = predictFromRF(zTrace, toeplitz, ...
    spatTempRF)
%PREDICTFROMRF   Predict response trace from stimulus and receptive field.

% INPUTS
% zTrce             [t], z-scored calcium trace sampled at stimulus 
%                   presentation times
% toeplitz          [t x pixels]; noise stimulus
% spatTempRF        [rows x columns (x sizes) x t_rf x ON/OFF], 2D masks of
%                   fitted ON and OFF subfields; in case of circle stimuli:
%                   masks for each circle size

% OUTPUTS
% prediction        [t], predicted trace
% explainedVar      double, explained variance

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