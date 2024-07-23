function writeKernelFitResults(fitResults, time, folder, type)
%WRITEKERNELFITRESULTS   Save results of kernel fitting procedure to
%explain calcium responses to visual stimuli

% INPUTS
% fitResults
%   .kernel         [k x ROIs], shape of kernel, normalized to max = 1
%   .amplitudes     [rep x stim x ROIs], kernel amplitudes for each repetition
%                   (rep) of each stimulus
%   .lags           [rep x stim x ROIs], kernel lag for each repetition of
%                   each stimulus
%   .prediction     [time x ROIs], predicted calcium trace
%   .pValue         [ROIs x 1], p-value of fit based on shift-test; set to 2
%                   if more than 5% of current set of shifted data
%                   resulted in better fit than non-shifted data
%   .R2             [ROIs x 1], explained variance of fit
% time
%   .kernel         [k x 1], time of kernel samples
%   .prediction     [time x 1], time of predicted trace samples
% folder            path to data of recording session
% type              str, stimulus type:
%                   'grating'/'gratingsDrifting'/'gratingsStatic'/'bars'

% time of kernel samples [t_k x 1]
writeNPY(time.kernel, fullfile(folder, ...
    sprintf('_ss_%sKernels.timestamps.npy', type)))
% shape of kernels [t_k x ROIs]
writeNPY(cat(2, fitResults.kernel), fullfile(folder, ...
    sprintf('_ss_%sKernels.dff.npy', type)))
% kernel amplitudes for each trial [repetition x stim x ROIs]
writeNPY(cat(3, fitResults.amplitudes), fullfile(folder, ...
    sprintf('_ss_%sTrials.amplitudes.npy', type)))
% p-value of fits [ROIs x 1]
writeNPY(cat(1, fitResults.pValue), fullfile(folder, ...
    sprintf('_ss_%sFits.pValue.npy', type)))
% explained variane of fits [ROIs x 1]
writeNPY(cat(1, fitResults.R2), fullfile(folder, ...
    sprintf('_ss_%sFits.R2.npy', type)))
% predicted calcium traces [t x ROIs]
writeNPY(cat(2, fitResults.prediction), fullfile(folder, ...
    sprintf('_ss_%sPredictions.dff.npy', type)))
% time of predicted trace samples [t x 1]
writeNPY(cat(2, time.prediction), fullfile(folder, ...
    sprintf('_ss_%sPredictions.timestamps.npy', type)))

if strcmp(type, 'bars')
    % optimal kernel shifts for each stimulus (for bars only) [repetition x
    % stim x ROIs]
    writeNPY(cat(3, fitResults.lags), fullfile(folder, ...
        sprintf('_ss_%sTrials.lags.npy', type)))
end