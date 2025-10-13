function [stimData_new, dirExp] = fixNonDirectionFeatures(stimData)

stimData_new = stimData;

nonDirFeatures = [stimData.contrasts, stimData.spatialFrequencies, ...
    stimData.temporalFrequencies];

% first, separate experiments, then find the one where only direction was 
% varied
ISIs = stimData.times(2:end, 1) - stimData.times(1:end-1, 2);
lastExpTrials = find(ISIs > 5);
if isempty(lastExpTrials)
    starts = 1;
    ends = size(stimData.times,1);
else
    starts = [1; lastExpTrials+1];
    ends = [lastExpTrials; size(stimData.times,1)];
end
dirExp = [];
stimIDs = [];
uniqueFeatures = [];
for k = 1:length(starts)
    ids = unique(stimData.ids(starts(k):ends(k)));
    uniqueStims = unique(nonDirFeatures(ids,:), "rows");
    if size(uniqueStims,1) == 1
        dirExp = [dirExp; k];
        stimIDs = [stimIDs; ids];
        uniqueFeatures = [uniqueFeatures; uniqueStims];
    end
end

if isempty(dirExp)
    fprintf('\n    NOTE: no experiment with correct stimulus features!\n\n')
elseif length(dirExp) > 1
    % if several direction experiments were performed, find the one where
    % other stimulus features are closest to the standard
    deviance = sum(abs((uniqueFeatures - [1 0.08 2]) ./ [1 0.08 2]), 2);
    [~, best] = min(deviance);
    if any(deviance == 0)
        dirExp = dirExp(deviance == 0);
    else
        dirExp = dirExp(best);
    end
end

% Update stimData_new to only include trials with most common features
indTrials = [];
for k = 1:length(dirExp)
    indTrials = [indTrials; (starts(dirExp(k)):ends(dirExp(k)))'];
end
stimData_new.times = stimData.times(indTrials,:);
stimData_new.ids = stimData.ids(indTrials);
rejected = find(~ismember((1:length(stimData.contrasts))', ...
    unique(stimData_new.ids)));
stimData_new.directions(rejected) = NaN;
stimData_new.contrasts(rejected) = NaN;
stimData_new.spatialFrequencies(rejected) = NaN;
stimData_new.temporalFrequencies(rejected) = NaN;