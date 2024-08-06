function mergedImage = mergeImageWithMasks(image, map, colors)
%MERGEIMAGEWITHMASKS   Generate one RGC image showing the mean frame with
%ROI masks marked.

% INPUTS
% image     [rows x columns], mean image
% map       [rows x columns], image with ROI masks, output from 
%           spatial.getROIMaskImage
% colors    [ROIs x 3], colors of ROI masks, output from
%           spatial.plotROIMaskImage

% OUTPUTS
% mergedImage [rows x columns x 3], RGB image showing input image and
%           transparent ROI masks plotted on top

% map ROI masks into RGB image
mapRGB = zeros([size(map), 3]);
mapRGB = reshape(mapRGB, numel(map), 3);
for col = 1:max(map(:))
    ind = find(map == col);
    mapRGB(ind,:) = repmat(colors(col+1,:), length(ind), 1);
end

% normalize image to range [0 1]
image = image - min(image(:));
image = image ./ max(image(:));
% transform image into RGC space
mergedImage = repmat(image, 1, 1, 3);
% add masks onto image
mergedImage = reshape(mergedImage, numel(image), 3);
mergedImage = mergedImage + 0.2 * mapRGB;
colored = map > 0;
colored = colored(:);
mergedImage(colored,:) = mergedImage(colored,:) ./ 1.2;
mergedImage = reshape(mergedImage, [size(image) 3]);