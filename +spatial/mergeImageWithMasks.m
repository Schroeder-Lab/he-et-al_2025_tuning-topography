function mergedImage = mergeImageWithMasks(image, map, colors)
% image     mean image with size of x- and y-range, all entries in the
%           range [0 1]
% mask      output from spatial.plotCellMasks
% colors    output from spatial.plotCellMasks

mapRGB = zeros([size(map), 3]);
mapRGB = reshape(mapRGB, numel(map), 3);
for col = 1:max(map(:))
    ind = find(map == col);
    mapRGB(ind,:) = repmat(colors(col+1,:), length(ind), 1);
end

image = image - min(image(:));
image = image ./ max(image(:));
mergedImage = repmat(image, 1, 1, 3);
mergedImage = reshape(mergedImage, numel(image), 3);
mergedImage = mergedImage + 0.2 * mapRGB;
colored = map > 0;
colored = colored(:);
mergedImage(colored,:) = mergedImage(colored,:) ./ 1.2;
mergedImage = reshape(mergedImage, [size(image) 3]);