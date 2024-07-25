function fig = plotPrefDiffVsDist(distances, prefDiffs, prefDiffsPermuted, ...
    binSize, stepSize)

fig = [];
minSamples = 6;

numPerm = size(prefDiffsPermuted,2);

binCentres = 0:stepSize:max(distances);
binValid = true(length(binCentres),1);
meanDist = NaN(length(binCentres),1);
meanDiff = NaN(length(binCentres),1);
semDiff = NaN(length(binCentres),1);
meanDiffPerm = NaN(length(binCentres), numPerm);
for bin = 1:length(binCentres)
    ind = distances > binCentres(bin)-binSize/2 & ...
        distances <= binCentres(bin)+binSize/2;
    meanDist(bin) = mean(distances(ind));
    if sum(ind) < minSamples
        binValid(bin) = false;
    end

    meanDiff(bin) = mean(prefDiffs(ind), 'omitnan');
    semDiff(bin) = std(prefDiffs(ind), 'omitnan') ./ sqrt(sum(ind));
    meanDiffPerm(bin,: ) = mean(prefDiffsPermuted(ind,:), 1, 'omitnan');
end
if sum(binValid) < 2
    return
end
meanDist = meanDist(binValid);
meanDiff = meanDiff(binValid);
meanDiffNull = prctile(meanDiffPerm(binValid,:), [2.5 50 97.5], 2);
x = ceil(meanDist(1):2:floor(meanDist(end)));
y = interp1(meanDist, meanDiff, x);
ySm = smoothdata(y, 'gaussian', 7);
yNull = interp1(meanDist, meanDiffNull, x);
yNullSm = smoothdata(yNull, 1, 'gaussian', 7);

h = [0 0];
fig = figure;
plot(distances, prefDiffs, 'k.')
hold on
h(1) = plot(x, ySm, 'r-', 'LineWidth', 2);
fill([x flip(x)], [yNullSm(:,1); flip(yNullSm(:,3))], 'b', ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none')
h(2) = plot(x, yNullSm(:,2), 'b-', 'LineWidth', 2);
legend(h, 'mean', 'mean of null data')
xlabel('Distance (um)')
ylabel('\DeltaPreference')