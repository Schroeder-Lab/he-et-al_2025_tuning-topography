function [ds_trans, os_long, os_lat] = getDsOsAxes(pitchAngle)

% Email from Shai Sabbah, 03.04.2025, all axes in [azimuth, elevation], 
% all axes provided for right eye!
% Azimuth: + -> towards right
% Elevation: + -> upwards
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

% NOTE: 
% For all vectors, we are interested in the left eye (mirror along saggital
% plane: multiply with [-1 1]).
% Vectors for direction coding point to optimal mouse-centred heading
% direction (see Extended Data Figure 9e in Sabbah et al (2017)). We,
% however, want the flow direction, which points into the opposite
% direction (add [180, 0], multiply with [1 -1]).

% Translation vectors for DSGCs from Sabbah et al (2017), here for left eye
% [azimuth, elevation]
ds_N = [180 -15]; % advance cells (optic flow), opposite to [-0 15] (black)
ds_T = [0 11]; % retreat cells, opposite to [-180 -11] (blue)
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

if nargin > 0
    [ds_trans(:,1), ds_trans(:,2)] = ...
        algebra.rotateY(ds_trans(:,1), ds_trans(:,2), pitchAngle);
    [os_long(:,1), os_long(:,2)] = ...
        algebra.rotateY(os_long(:,1), os_long(:,2), pitchAngle);
    [os_lat(:,1), os_lat(:,2)] = ...
        algebra.rotateY(os_lat(:,1), os_lat(:,2), pitchAngle);
end