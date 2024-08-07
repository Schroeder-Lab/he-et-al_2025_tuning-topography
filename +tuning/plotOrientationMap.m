function fig = plotOrientationMap(angles, isTuned, type, masks, fovPix, fovMicron)
%PLOTORIENTATIONMAP   Plot ROI masks within imaging FOV with preferred
%direction/orientation color coded.

% INPUTS
% angles    [ROIs x 1], preferred directions or orientations of ROIs, NaN
%           if not responsive (kernel fit not good)
% isTuned   [ROIs x 1], true if significantly tuned, otherwise false
% type      'ori' or 'dir'
% masks     [ROIs x pixels], location of ROI masks indexing into foxPix
% fovPix    [1 x 2]; size of frame in pixels (vertical, horizontal)
% fovMicron [1 x 2], size of frame in microns (vertical, horizontal)

% OUTPUTS
% fig       int, handle to figure

% colors used for plotting: white: background, light gray: not responsive,
% darker gray: not tuned, colors: represent angle
colors = [hsv(360); 0.6 0.6 0.6; 0.8 0.8 0.8; 1 1 1];
% create image of FOV (all black)
map = ones(fovPix) .* 363;
% set directions/orientations to range 1-360/1-180
if strcmp(type, 'ori')
    angles = round(mod(angles(:,1),180)*2);
    str = 'Orientation';
    labels = 0:15:180;
elseif strcmp(type, 'dir')
    angles = round(mod(angles(:,1),360));
    str = 'Direction';
    labels = 0:30:360;
end
angles(angles == 0) = 360;
% set untuned ROIs to light gray
angles(~isTuned & ~isnan(angles)) = 361;
% set unresponsive ROIs to dark gray
angles(isnan(angles)) = 362;
% for each ROI, mark mask in color of preferred direction/orientation
for roi = 1:size(angles,1)
    m = masks(roi,:);
    m(isnan(m)) = [];
    map(m) = angles(roi);
end

fig = figure('Position', [525 150 945 780]);
% plot masks
ax(1) = subplot(1,2,1);
image([0 fovMicron(2)], [0 fovMicron(1)], map)
colormap(colors);
axis image
xlabel('FOV position (um)')
title([str ' map'])
% plot color map
ax(2) = subplot(1,2,2);
image((1:360)')
colormap(colors);
ylabel(str)
% rescale plots
set(ax(2), 'Position', [0.85 0.11 0.05 0.82], 'YAxisLocation', 'right', ...
    'YTick', 1:30:360, 'YTickLabel', labels, 'XTick', [], 'box', 'off')
set(ax(1), 'Position', [0.13 0.11 0.7 0.82])
axis(ax(1), 'equal')