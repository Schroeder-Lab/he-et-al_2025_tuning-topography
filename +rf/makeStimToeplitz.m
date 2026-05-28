function st = makeStimToeplitz(stim, rfBins)
%MAKESTIMTOEPLITZ   Generate toeplitz matrix for visual noise stimuli.

% INPUTS
%   stim                [time x pixels]; noise stimulus
%   rfBins              [1 x RFframes]; time of spatio-temporal RF relative
%                       to response in number of stimulus frames

% OUTPUTS
%   stim                [time x pixels]; toeplitz matrix

% generate toeplitz matrix for stimuli: [time x pixels]
% each row holds all pixels at current and previous time points:
% [[all pixels at t=0], [all pixels at t=-1], ...]
% each column is time series of that particular pixel

% concatinate time shifted stimulus blocks; for each time point there
% is a stimulus block for lag=0, another for lag=-1, another for lag=-2,...
st = [];
for t = 1:length(rfBins)
    st = [st, ...
        [zeros(max(0, rfBins(1) - 1 + t), size(stim, 2)); ...
        stim(max(1, 2 - rfBins(1) - t) : end - rfBins(1) - t + 1, :)]];
end