function [direction, orientation] = determineTuning(amplitudes, ...
    stimDirections, numShuffles)
%DETERMINETUNING   Determine direction and orientation preference,
%selectivity and significance of one unit.

% INPUTS
% amplitudes        [repetitions x stim], response amplitudes to multiple
%                   repetitions of the same stimuli
% stimDirections    [stim x 1], direction of stimuli in degrees
% numShuffles       int, number of permutations for significance 
%                   (permutation) test

% OUTPUTS
% direction
%   .preference     double, preferred direction (direction of vector sum)
%   .selectivity    double, direction selectivity (length of vector sum)
%   .responseSign   1 or -1, sign of mean response amplitude (across all
%                   stimuli)
%   .pValue         p-value of permutation test
% orientation
%   .preference     double, preferred orientation (half direction of vector sum)
%   .selectivity    double, orientation selectivity (length of vector sum)
%   .responseSign   1 or -1, sign of mean response amplitude (across all
%                   stimuli)
%   .pValue         p-value of permutation test

% direction tuning
[pref, sel, sgn] = tuning.vectorAveraging(amplitudes, stimDirections);
direction.preference = pref;
direction.selectivity = sel;
direction.responseSign = sgn;

% orientation tuning
% treat opposite movement directions as the same orientation
stimOrientations = mod(stimDirections, 180);
% determine preference, selectivity, significance
[pref, sel, sgn] = tuning.vectorAveraging(amplitudes, 2*stimOrientations);
orientation.preference = pref/2;
orientation.selectivity = sel;
orientation.responseSign = sgn;

% do permutation test
numTrials = numel(amplitudes);
numRep = size(amplitudes,1);
dirSelectivities = NaN(numShuffles,1);
oriSelectivities = NaN(numShuffles,1);
for k = 1:numShuffles
    indShuffled = randperm(numTrials);
    ampsShuffled = reshape(amplitudes(indShuffled), numRep, []);
    [~, dirSelectivities(k)] = tuning.vectorAveraging(ampsShuffled, ...
        stimDirections);
    [~, oriSelectivities(k)] = tuning.vectorAveraging(ampsShuffled, ...
        2*stimOrientations);
end

% determine p-values based on permuted data
direction.pValue = sum(dirSelectivities > direction.selectivity) / numShuffles;
orientation.pValue = sum(oriSelectivities > orientation.selectivity) / numShuffles;