function data = getStimResponseFits(folder, type)
%GETSTIMULUSRESPONSEFITS  Load results of kernel fitting procedure to
%explain calcium responses to visual stimuli

% INPUTS
% folder            path to data of recording session
% type              str, stimulus type:
%                   'grating'/'gratingsDrifting'/'gratingsStatic'/'bars'

% OUTPUTS
% data
%   .time_kernel    [k x 1], time of kernel samples
%   .kernel         [k x ROIs], shape of kernel, normalized to max = 1
%   .amplitudes     [rep x stim x ROIs], kernel amplitudes for each repetition
%                   (rep) of each stimulus
%   .time_prediction [time x 1], time of predicted trace samples
%   .prediction     [time x ROIs], predicted calcium trace
%   .pValue         [ROIs x 1], p-value of fit based on shift-test; set to 2
%                   if more than 5% of current set of shifted data
%                   resulted in better fit than non-shifted data
%   .R2             [ROIs x 1], explained variance of fit
%   (.lags)         (only for bars)
%                   [rep x stim x ROIs], kernel lag for each repetition of
%                   each stimulus

data.time_kernel = readNPY(fullfile(folder, sprintf('_ss_%sKernels.timestamps.npy', type)));
data.kernel = readNPY(fullfile(folder, sprintf('_ss_%sKernels.dff.npy', type)));
data.amplitudes = readNPY(fullfile(folder, sprintf('_ss_%sTrials.amplitudes.npy', type)));
data.time_prediction = readNPY(fullfile(folder, sprintf('_ss_%sPredictions.timestamps.npy', type)));
data.prediction = readNPY(fullfile(folder, sprintf('_ss_%sPredictions.dff.npy', type)));
data.pValue = readNPY(fullfile(folder, sprintf('_ss_%sFits.pValue.npy', type)));
data.R2 = readNPY(fullfile(folder, sprintf('_ss_%sFits.R2.npy', type)));

if strcmp(type, 'bars')
    data.lags = readNPY(fullfile(folder, sprintf('_ss_%sTrials.lags.npy', type)));
end