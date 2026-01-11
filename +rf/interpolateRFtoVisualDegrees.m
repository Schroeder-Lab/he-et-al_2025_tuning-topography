function [rf_visDeg, x2, y2] = interpolateRFtoVisualDegrees(rf, stimPos)
%INTERPOLATERFTOVISUALDEGREES   Sample RF profile so that grid points are
%separated by one visual degree.

% INPUTS
% rf        [rows x columns], 2D RF profile
% stimPos   [left right top bottom], limits of stimulus in visual degrees

% OUTPUTS
% rf_visDeg [newRows x newColumns], 2D RF profile resampled
% x2        [newRows x newColumns], azimuth of each RF sample point
% y2        [newRows x newColumns], elevation of each RF sample point

[x0, y0] = meshgrid((1:size(rf,2))-0.5, (1:size(rf,1))-0.5);
% [x0, y0] = meshgrid(1:size(rf,2), 1:size(rf,1));

% vectors x1 and y1 specify gridpoints with distance of 1
% degree (diff(stimPos(...))); values match position of
% gridpoints in pixels of stimulus row/column
x1 = linspace(0, size(rf,2), diff(stimPos(1:2)));
y1 = linspace(0, size(rf,1), -diff(stimPos(3:4)));
% x1 = linspace(0.5, size(rf,2)+0.5, diff(stimPos(1:2)));
% y1 = linspace(0.5, size(rf,1)+0.5, -diff(stimPos(3:4)));

% need to delete pixel values outside given stimulus pixels,
% so we can use interpolation (rather than extrapolation) when
% mapping the RF from stimulus pixels to visual degrees
x2 = x1;
x2(x1<x0(1) | x1>x0(end)) = [];
% x2(x1<1 | x1>size(rf,2)) = [];
y2 = y1;
y2(y1<y0(1) | y1>y0(end)) = [];
% y2(y1<1 | y1>size(rf,1)) = [];
[x2, y2] = meshgrid(x2, y2);

% interpolate given RF to new grid
rf_visDeg = interp2(x0, y0, rf, x2, y2);