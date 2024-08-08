function [prefDiffsRelative, meanDist] = getPrefDiffsRelativeNull(distances, ...
    prefDiffs, prefDiffsPermuted, binSize, stepSize)
%GETPREFDIFFSRELATIVENULL   Determine mean preference difference in pairs
%depending on their distance; measure differences relative to null
%distribution of differences.

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
% prefDiffsRelative [bins], preference differences per bin, relative to
%                   null distribution (in range [0 1])
% meanDist          [bins], mean distance across pairs per bin

[meanDist, meanDiff, meanDiffPerm] = spatial.getMeanPrefDiffPerDistBin( ...
    distances, prefDiffs, prefDiffsPermuted, binSize, stepSize);

% express actual preference differences in terms of percentiles in null
% distribution
if all(isnan(prefDiffs))
    prefDiffsRelative = NaN(size(meanDiff));
else
    prefDiffsRelative = sum(meanDiff > meanDiffPerm, 2) ./ ...
        size(meanDiffPerm,2);
end