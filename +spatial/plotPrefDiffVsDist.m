function fig = plotPrefDiffVsDist(distances, prefDiffs, prefDiffsPermuted, ...
    binSize, stepSize, doSmooth)
%PLOTPREFDIFFVSDIST   Plot brain distances versus differences in
%direction/orientation preferences.

% INPUTS
% distance      [pairs x 1], pairwise distances
% prefDiffs     [pairs x 1], pairwise tuning differences
% prefDiffsPermuted [pairs x permutations], tuning differences of permuted
%               ROIs
% binSize       int, bin size to determine average differences
% stepSize      int, step size for moving average

% OUTPUTS
% fig       int, handle to figure

fig = [];

[meanDist, meanDiff, meanDiffPerm] = spatial.getMeanPrefDiffPerDistBin( ...
    distances, prefDiffs, prefDiffsPermuted, binSize, stepSize);
if all(isnan(meanDist))
    return
end
% determine median and confidence interval of null distribution (permuted
% data)
meanDiffNull = prctile(meanDiffPerm, [2.5 50 97.5], 2);
% interpolate data to equal-distant brain distances
x = ceil(meanDist(1):2:floor(meanDist(end)));
y = interp1(meanDist, meanDiff, x);
yNull = interp1(meanDist, meanDiffNull, x);
if doSmooth
    % smooth moving averages
    ySm = smoothdata(y, 'gaussian', 7);
    yNullSm = smoothdata(yNull, 1, 'gaussian', 7);
else
    ySm = y;
    yNullSm = yNull;
end

h = [0 0];
fig = figure;
% scatterplot
plot(distances, prefDiffs, 'k.')
hold on
% moving average of real data
h(1) = plot(x, ySm, 'r-', 'LineWidth', 2);
% moving average of null data and confidence interval
fill([x flip(x)], [yNullSm(:,1); flip(yNullSm(:,3))], 'b', ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none')
h(2) = plot(x, yNullSm(:,2), 'b-', 'LineWidth', 2);
legend(h, 'mean', 'mean of null data')
xlabel('Distance (um)')
ylabel('\DeltaPreference')