function writeKernelFitResults(fitResults, time, folder, type)

if nargin < 2
    type = '';
end

switch type
    case 'drifting'
        str = 'gratingsDrifting';
    case 'static'
        str = 'gratingsStatic';
    otherwise
        str = 'grating';
end

writeNPY(time.kernel, fullfile(folder, ...
    sprintf('_ss_%sKernels.timestamps.npy', str)))
writeNPY(cat(2, fitResults.kernel), fullfile(folder, ...
    sprintf('_ss_%sKernels.dff.npy', str)))
writeNPY(cat(3, fitResults.amplitudes), fullfile(folder, ...
    sprintf('_ss_%sTrials.amplitudes.npy', str)))
writeNPY(cat(1, fitResults.pValue), fullfile(folder, ...
    sprintf('_ss_%sFits.pValue.npy', str)))
writeNPY(cat(1, fitResults.R2), fullfile(folder, ...
    sprintf('_ss_%sFits.R2.npy', str)))
writeNPY(cat(2, fitResults.prediction), fullfile(folder, ...
    sprintf('_ss_%sPredictions.dff.npy', str)))
writeNPY(cat(2, time.prediction), fullfile(folder, ...
    sprintf('_ss_%sPredictions.timestamps.npy', str)))