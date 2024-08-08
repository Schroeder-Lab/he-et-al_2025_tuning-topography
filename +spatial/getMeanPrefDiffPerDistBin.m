function [meanDist, meanDiff, meanDiffPerm] = ...
    getMeanPrefDiffPerDistBin(distances, prefDiffs, prefDiffsPermuted, ...
    binSize, stepSize)
%GETMEANPREFDIFFPERDISTBIN   Bin pairs of units according to their brain
%distance, then calculate mean difference between preferences across those
%pairs (and null distribution).

% INPUTS
% distances         [pairs], distances of pairs in microns
% prefDiffs         [pairs], difference between direction/orientation
%                   preferences for each pair
% prefDiffsPermuted [pairs x permutations], preference differences after
%                   permuting units
% binSize           double, range of brain distance to consider in each bin
%                   (in microns)
% stepSize          double, distance between bins (in microns)

% OUTPUTS
% meanDist          [bins], mean distance across pairs per bin
% meanDiff          [bins], mean preference difference across pairs per bin
% meanDiffPerm      [bins x permutations], mean preference difference
%                   across pairs per bin, calculated for each permutation

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
    meanDiff(bin) = mean(prefDiffs(ind), 'omitnan');
    meanDiffPerm(bin,:) = mean(prefDiffsPermuted(ind,:), 1, 'omitnan');
    
    % check that bin contains sufficient samples
    if sum(ind) < minSamples
        binValid(bin) = false;
    end
end
% ignore bins with insufficient pairs
meanDist = meanDist(binValid);
meanDiff = meanDiff(binValid);
meanDiffPerm = meanDiffPerm(binValid,:);