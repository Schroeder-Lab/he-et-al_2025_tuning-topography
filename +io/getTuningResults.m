function [dirTuning, oriTuning] = getTuningResults(folder, type)
%GETTUNINGDATA  Load direction and orientation preferences, selectivities
%and significance results.

% INPUTS
% folder            path to data of recording session
% type              str, stimulus type:
%                   'grating'/'gratingsDrifting'/'gratingsStatic'/'bars'

% OUTPUTS
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

if any(strcmp(type, {'gratingsDrifting', 'bars'}))
    dirTuning.preference = readNPY(fullfile(folder, ...
        sprintf('_ss_%sTuning.directionPreference.npy', type)));
    dirTuning.selectivity = readNPY(fullfile(folder, ...
        sprintf('_ss_%sTuning.directionSelectivity.npy', type)));
    dirTuning.pValue = readNPY(fullfile(folder, ...
        sprintf('_ss_%sTuning.directionPValue.npy', type)));
    dirTuning.responseSign = readNPY(fullfile(folder, ...
        sprintf('_ss_%sTuning.directionSign.npy', type)));
end
oriTuning.preference = readNPY(fullfile(folder, ...
    sprintf('_ss_%sTuning.orientationPreference.npy', type)));
oriTuning.selectivity = readNPY(fullfile(folder, ...
    sprintf('_ss_%sTuning.orientationSelectivity.npy', type)));
oriTuning.pValue = readNPY(fullfile(folder, ...
    sprintf('_ss_%sTuning.orientationPValue.npy', type)));
oriTuning.responseSign = readNPY(fullfile(folder, ...
    sprintf('_ss_%sTuning.orientationSign.npy', type)));