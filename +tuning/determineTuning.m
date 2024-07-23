function [direction, orientation] = determineTuning(amplitudes, ...
    stimDirections, numShuffles)

% direction tuning
[pref, sel, sgn] = tuning.vectorAveraging(amplitudes, stimDirections);
direction.preference = pref;
direction.selectivity = sel;
direction.responseSign = sgn;

% orientation tuning
stimOrientations = mod(stimDirections, 180);
[pref, sel, sgn] = tuning.vectorAveraging(amplitudes, 2*stimOrientations);
orientation.preference = pref/2;
orientation.selectivity = sel;
orientation.responseSign = sgn;

% do shuffle test
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

% determine p-values based on shuffled data
direction.pValue = sum(dirSelectivities > direction.selectivity) / numShuffles;
orientation.pValue = sum(oriSelectivities > orientation.selectivity) / numShuffles;