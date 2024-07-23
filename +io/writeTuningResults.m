function writeTuningResults(dirTuning, oriTuning, folder, type)
%WRITETUNINGRESULTS   Save results of direction and orientation tuning

% INPUTS
% dirTuning
%   .preference     [ROIs x 1], preferred directions
%   .selectivity    [ROIs x 1], direction selectivities
%   .pValue         [ROIs x 1], p-values for selectivities based on permutation-test
%   .responseSign   [ROIs x 1], signs of average response (per ROI)
% oriTuning
%   .preference     [ROIs x 1], preferred orientations
%   .selectivity    [ROIs x 1], orientation selectivities
%   .pValue         [ROIs x 1], p-values for selectivities based on permutation-test
%   .responseSign   [ROIs x 1], signs of average response (per ROI)
% folder            path to data of recording session
% type              str, stimulus type:
%                   'grating'/'gratingsDrifting'/'gratingsStatic'/'bars'

if any(strcmp(type, {'gratingsDrifting', 'bars'}))
    writeNPY(cat(1, dirTuning.preference), fullfile(folder, ...
        sprintf('_ss_%sTuning.directionPreference.npy', type)))
    writeNPY(cat(1, dirTuning.selectivity), fullfile(folder, ...
        sprintf('_ss_%sTuning.directionSelectivity.npy', type)))
    writeNPY(cat(1, dirTuning.pValue), fullfile(folder, ...
        sprintf('_ss_%sTuning.directionPValue.npy', type)))
    writeNPY(cat(1, dirTuning.responseSign), fullfile(folder, ...
        sprintf('_ss_%sTuning.directionSign.npy', type)))
end
writeNPY(cat(1, oriTuning.preference), fullfile(folder, ...
    sprintf('_ss_%sTuning.orientationPreference.npy', type)))
writeNPY(cat(1, oriTuning.selectivity), fullfile(folder, ...
    sprintf('_ss_%sTuning.orientationSelectivity.npy', type)))
writeNPY(cat(1, oriTuning.pValue), fullfile(folder, ...
    sprintf('_ss_%sTuning.orientationPValue.npy', type)))
writeNPY(cat(1, oriTuning.responseSign), fullfile(folder, ...
    sprintf('_ss_%sTuning.orientationSign.npy', type)))