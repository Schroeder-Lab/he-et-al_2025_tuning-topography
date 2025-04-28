function [ds_trans, os_long, os_lat] = getDsOsAxes()

% Email from Shai Sabbah, 03.04.2025, all axes in [azimuth, elevation], 
% all axes provided for right eye!
% For DSGCs:
% N, 0, 15
% T, 180, -11
% D, -113, 83
% V, 74, -88
% For OSGCs:
% Longitudinal system:
% H, -7, 1
% V, 82, -64
% Latitudinal system:
% H, 105, -68
% V, -53, -6 (opposite: [127 6])

% Translation vectors for DSGCs from Sabbah et al (2017), here for left eye
% [azimuth, elevation]
ds_N = [180 -15]; % advance cells, opposite to [0 15] (black)
ds_T = [0 11]; % retreat cells, opposite to [180 -11] (blue)
ds_D = [293 -83]; % rise cells, opposite to [113 83] (red)
ds_V = [106 88]; % fall cells, opposite to [-74 -88] (magenta)
ds_trans = [ds_N; ds_T; ds_D; ds_V];

% Axes for OSGCs from Laniado et al (2025), here for left eye
% longitudinal cells
os_long_H = [7 1]; % black
os_long_V = [-82 -64]; % blue
os_long = [os_long_H; os_long_V];
% latitudinal cells
os_lat_H = [-105 -68]; % red
os_lat_V = [53 -6]; % magenta (opposite: [-127 6])
os_lat = [os_lat_H; os_lat_V];