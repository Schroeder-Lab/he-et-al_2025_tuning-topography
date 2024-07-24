function fig = plotOrientationMap(angles, isTuned, type, masks, fovPix, fovMicron)

% neuronOrientations    [neurons x 2]; 1st col: pref. orientation or
%                       direction, 2nd col.: goodness of tuning, e.g. adj.
%                       R-Squared
% type                  'ori' or 'dir'
% ROImaps               {neurons x 1}; each entry contains pixels of ROI
% mapsSize              [1 x 2]; size of frame ([validY validX])
% minWeight             minimum goodness of tuning to show orientation,
%                       otherwise ROI will be gray

% map = zeros([fovPix,3]);
map = ones(fovPix) .* 363;
colors = [hsv(360); 0.4 0.4 0.4; 0.2 0.2 0.2; 0 0 0];
if strcmp(type, 'ori')
    angles = round(mod(angles(:,1),180)*2);
    str = 'Orientation';
elseif strcmp(type, 'dir')
    angles = round(mod(angles(:,1),360));
    str = 'Direction';
end
angles(angles == 0) = 360;
angles(~isTuned & ~isnan(angles)) = 361;
angles(isnan(angles)) = 362;
for roi = 1:size(angles,1)
    m = masks(roi,:);
    m(isnan(m)) = [];
    % map(m) = c(roi,1);
    % map(prod(fovPix) + m) = c(roi,2);
    % map(2*prod(fovPix) + m) = c(roi,3);
    map(m) = angles(roi);
end

fig = figure('Position', [525 150 945 780]);
ax(1) = subplot(1,2,1);
% image(map)
image([0 fovMicron(2)], [0 fovMicron(1)], map)
colormap(colors);
axis image
xlabel('FOV position (um)')
title([str ' map'])
ax(2) = subplot(1,2,2);
image((1:360)')
colormap(colors);
ylabel(str)
set(ax(2), 'Position', [0.85 0.11 0.05 0.82], 'YAxisLocation', 'right', ...
    'YTick', 1:30:360, 'YTickLabel', 0:15:179, 'XTick', [], 'box', 'off')
set(ax(1), 'Position', [0.13 0.11 0.7 0.82])
axis(ax(1), 'equal')