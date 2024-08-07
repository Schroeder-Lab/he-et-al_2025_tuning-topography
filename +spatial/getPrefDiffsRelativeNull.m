function [prefDiffsRelative, meanDist] = getPrefDiffsRelativeNull(distances, ...
    prefDiffs, prefDiffsPermuted, binSize, stepSize)

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