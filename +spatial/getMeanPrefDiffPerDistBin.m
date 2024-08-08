function [meanDist, meanDiff, meanDiffPerm] = ...
    getMeanPrefDiffPerDistBin(distances, prefDiffs, prefDiffsPermuted, ...
    binSize, stepSize)

minSamples = 6; % minimum number of pairs per bin for moving average
numPerm = size(prefDiffsPermuted,2);
binEdges = 0:stepSize:max(distances);
binValid = true(length(binEdges),1);
meanDist = NaN(length(binEdges),1);
meanDiff = NaN(length(binEdges),1);
meanDiffPerm = NaN(length(binEdges), numPerm);
% loop over bins to determine moving average (of tuning differences and
% brain distances)
for bin = 1:length(binEdges)
    ind = distances > binEdges(bin)-binSize/2 & ...
        distances <= binEdges(bin)+binSize/2;
    meanDist(bin) = mean(distances(ind));
    % check that bin contains sufficient samples
    if sum(ind) < minSamples
        binValid(bin) = false;
    end

    meanDiff(bin) = mean(prefDiffs(ind), 'omitnan');
    meanDiffPerm(bin,:) = mean(prefDiffsPermuted(ind,:), 1, 'omitnan');
end
% ignore bins with insufficient pairs
meanDist = meanDist(binValid);
meanDiff = meanDiff(binValid);
meanDiffPerm = meanDiffPerm(binValid,:);