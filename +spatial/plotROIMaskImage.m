function colors = plotROIMaskImage(map, masks, plotIDs)

numROIs = size(masks,1);
colors = [0 0 0; jet(numROIs)];

figure('Position', [680 50 1050 945])
imagesc(map)
colormap(colors)
axis image

if plotIDs
    centroids = regionprops(map, "Centroid");
    indWhite = (1:numROIs)';
    indWhite = indWhite < 0.28 * numROIs | indWhite > 0.65 * numROIs;
    valid = ~all(isnan(masks), 2);
    c = cat(1, centroids(indWhite & valid).Centroid);
    if ~isempty(c)
        text(c(:,1), c(:,2), num2str(find(indWhite & valid)), ...
            'Color', [1 1 1], ...
            'HorizontalAlignment', 'center', 'FontWeight', 'bold')
    end
    c = cat(1, centroids(~indWhite & valid).Centroid);
    if ~isempty(c)
        text(c(:,1), c(:,2), num2str(find(~indWhite & valid)), ...
            'Color', [1 1 1].*0.3, ...
            'HorizontalAlignment', 'center', 'FontWeight', 'bold')
    end
end