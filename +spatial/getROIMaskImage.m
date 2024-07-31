function map = getROIMaskImage(masks, fovPix, fovBoundaries)

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