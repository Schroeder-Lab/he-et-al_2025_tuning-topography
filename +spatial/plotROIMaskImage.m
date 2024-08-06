function colors = plotROIMaskImage(map, masks, plotIDs)
%PLOTROIMASKIMAGE   Plot patches representing ROI masks with/without ROI
%ID.

% INPUTS
% map       [rows x columns], image with ROI masks, output from 
%           spatial.getROIMaskImage
% masks     [ROIs x pixels], linear indices of each ROI mask within image
%           of same size as map
% plotIDs   true/false, if true add ROI IDs as text to masks

% OUTPUTS
% colors    colormap used to plot map

numROIs = size(masks,1);
colors = [0 0 0; jet(numROIs)];

figure('Position', [680 50 1050 945])
imagesc(map)
colormap(colors)
axis image

if plotIDs
    % get centres of ROI masks
    centroids = regionprops(map, "Centroid");
    % get text colour (white or dark grey)
    indWhite = (1:numROIs)';
    indWhite = indWhite < 0.28 * numROIs | indWhite > 0.65 * numROIs;
    % plot IDs
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