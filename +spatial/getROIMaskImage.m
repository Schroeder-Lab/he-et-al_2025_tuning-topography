function map = getROIMaskImage(masks, fovPix, fovBoundaries)
%GETROIMASKIMAGE   Create image with ROI masks from mask indices.

% INPUTS
% masks         [ROIs x pixels], linear indices of each ROI mask within image
%               of same size as map
% fovPix        [1 x 2], (numRows, numColumns) of field of view to place
%               ROI masks
% fovBoundaries [1 x 4], (top bottom left right)
%               pixels relative to full imaged FOV, pixels outside these
%               boundaries were disregarded for ROI detection as the edges
%               were not always visible

% OUTPUTS
% map           [rows x columns], image with ROI masks

map = zeros(fovPix);
% for each ROI, mark mask pixels with ROI ID (row in masks)
for roi = 1:size(masks,1)
    m = masks(roi,:);
    m(isnan(m)) = [];
    map(m) = roi;
end

if diff(fovBoundaries([1 2]))+1 ~= fovPix(1) || ...
        diff(fovBoundaries([3 4]))+1 ~= fovPix(2)
    map = map(fovBoundaries(1):fovBoundaries(2), ...
        fovBoundaries(3):fovBoundaries(4));
end