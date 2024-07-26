function fig = plotPrefDiffVsDist(distances, prefDiffs, prefDiffsPermuted, ...
    binSize, stepSize)
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
minSamples = 6; % minimum number of pairs per bin for moving average

numPerm = size(prefDiffsPermuted,2);
binCentres = 0:stepSize:max(distances);
binValid = true(length(binCentres),1);
meanDist = NaN(length(binCentres),1);
meanDiff = NaN(length(binCentres),1);
meanDiffPerm = NaN(length(binCentres), numPerm);
% loop over bins to determine moving average (of tuning differences and
% brain distances)
for bin = 1:length(binCentres)
    ind = distances > binCentres(bin)-binSize/2 & ...
        distances <= binCentres(bin)+binSize/2;
    meanDist(bin) = mean(distances(ind));
    % check that bin contains sufficient samples
    if sum(ind) < minSamples
        binValid(bin) = false;
    end

    meanDiff(bin) = mean(prefDiffs(ind), 'omitnan');
    meanDiffPerm(bin,: ) = mean(prefDiffsPermuted(ind,:), 1, 'omitnan');
end
% if not enough data, return
if sum(binValid) < 2
    return
end
% ignore bins with insufficient pairs
meanDist = meanDist(binValid);
meanDiff = meanDiff(binValid);
% determine median and confidence interval of null distribution (permuted
% data)
meanDiffNull = prctile(meanDiffPerm(binValid,:), [2.5 50 97.5], 2);
% interpolate data to equal-distant brain distances
x = ceil(meanDist(1):2:floor(meanDist(end)));
y = interp1(meanDist, meanDiff, x);
yNull = interp1(meanDist, meanDiffNull, x);
% smooth moving averages
ySm = smoothdata(y, 'gaussian', 7);
yNullSm = smoothdata(yNull, 1, 'gaussian', 7);

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