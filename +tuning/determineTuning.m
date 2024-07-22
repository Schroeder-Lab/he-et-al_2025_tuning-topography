function [direction, orientation] = determineTuning(amplitudes, ...
    stimDirections, numShuffles)

% direction tuning
[pref, sel] = tuning.vectorAveraging(amplitudes, stimDirections);
direction.preference = pref;
direction.selectivity = sel;

% orientation tuning
stimOrientations = mod(stimDirections, 180);
[pref, sel] = tuning.vectorAveraging(amplitudes, stimOrientations);
orientation.preference = pref;
orientation.selectivity = sel;

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
        stimOrientations);
end

% determine p-values based on shuffled data
direction.pValue = sum(dirSelectivities > direction.selectivity) / numShuffles;
orientation.pValue = sum(oriSelectivities > orientation.selectivity) / numShuffles;