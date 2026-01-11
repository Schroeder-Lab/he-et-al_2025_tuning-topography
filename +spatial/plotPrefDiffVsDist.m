function fig = plotPrefDiffVsDist(distances, prefDiffs, prefDiffsPermuted, ...
    binSize, stepSize, doSmooth)
%PLOTPREFDIFFVSDIST   Plot brain distances versus differences in
%direction/orientation preferences.

% INPUTS
% distances     [pairs x 1], pairwise distances
% prefDiffs     [pairs x 1], pairwise tuning differences
% prefDiffsPermuted [pairs x permutations], tuning differences of permuted
%               ROIs
% binSize       int, bin size to determine average differences
% stepSize      int, step size for moving average
% doSmooth      logical, if true smooth the mean preference differences
%               across x-axis (brain distance)

% OUTPUTS
% fig           int, handle to figure

fig = [];

[meanDist, meanDiff, meanDiffPerm] = spatial.getMeanPrefDiffPerDistBin( ...
    distances, prefDiffs, prefDiffsPermuted, binSize, stepSize);
if length(meanDist)<2 || all(isnan(meanDist))
    return
end
% remove duplicates
[meanDist, ind] = unique(meanDist, 'stable');
meanDiff = meanDiff(ind);
meanDiffPerm = meanDiffPerm(ind,:);
% determine median and confidence interval of null distribution (permuted
% data)
meanDiffNull = prctile(meanDiffPerm, [2.5 50 97.5], 2);
% interpolate data to equal-distant brain distances
d = median(diff(meanDist));
if d < 0.5
    d = 0.1;
elseif d < 1.5
    d = 0.5;
else
    d = 2;
end
x = (ceil(meanDist(1)/d):floor(meanDist(end)/d))  .* d;
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
% density of data points
dens = ksdensity([distances, prefDiffs], [distances, prefDiffs]);
% scatterplot
scatter(distances, prefDiffs, 10, dens, "filled")
colormap(gca, flip(gray))
c = colorbar;
c.Label.String = 'Density';
hold on
% moving average of real data
h(1) = plot(x, ySm, 'r-', 'LineWidth', 2);
% moving average of null data and confidence interval
fill([x flip(x)], [yNullSm(:,1); flip(yNullSm(:,3))], 'c', ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none')
h(2) = plot(x, yNullSm(:,2), 'c-', 'LineWidth', 2);
l = legend(h, 'mean', 'mean of null data');
l.Box = "off";
set(gca, 'box', 'off')
xlabel('Distance (um)')
ylabel('\DeltaPreference')