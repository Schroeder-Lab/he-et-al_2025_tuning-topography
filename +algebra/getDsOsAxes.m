function [ds_trans, os_long, os_lat] = getDsOsAxes(pitchAngle)

% Email from Shai Sabbah, 03.04.2025, all axes in [azimuth, elevation], 
% all axes provided for right eye!
% Azimuth: + -> towards right
% Elevation: + -> upwards
% (as described in Sabbah et al (2017) under "Transforming retinal 
% coordinates into global mouse-centred coordinates")
% For DSGCs:
% N, 0, 15
% T, 180, -11
% D, -113, 83
% V, 74, -88
% Coordinates for OSGCs in that email were incorrect.
 
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
ds_D = [-67 -83]; % rise cells, opposite to [113 83] (red)
ds_V = [106 88]; % fall cells, opposite to [-74 -88] (magenta)
ds_trans = [ds_N; ds_T; ds_D; ds_V];

% Email from Shai Sabbah, 31.10.2025:
% For each singularity, the first coordinate is the azimuth, the second 
% coordinate is the elevation.
% The azimuth coordinate is zero in front of the mouse and increases when 
% moving counterclockwise (when viewing the mouse from above). The 
% elevation coordinate is zero at the zenith and increasing toward the 
% nadir.
% Provided numbers for left eye:
% 1. Translation / longitudes
% H: [187, 90]
% V: [98, 26]
% 2. Rotation / latitudes
% H: [74, 22]
% V: [232, 85]

% For all vectors, let's transform the coordinates into the conventional
% format (Azimuth: 0 -> frontal, + -> towards right, Elevation: 0 -> 
% horizon, + -> upwards)
% Longitudes: H: [173, 0],  V: [-98, 63]
% Latitudes:  H: [-74, 68], V: [128, 5]

% Axes for OSGCs from Laniado et al (2025), here for left eye
% longitudinal cells
os_long_H = [173 0]; % black
os_long_V = [-98 63]; % blue
os_long = [os_long_H; os_long_V];
% latitudinal cells
os_lat_H = [-74 68]; % red
os_lat_V = [128 5]; % magenta
os_lat = [os_lat_H; os_lat_V];

if nargin > 0
    [ds_trans(:,1), ds_trans(:,2)] = ...
        algebra.rotateY(ds_trans(:,1), ds_trans(:,2), pitchAngle);
    [os_long(:,1), os_long(:,2)] = ...
        algebra.rotateY(os_long(:,1), os_long(:,2), pitchAngle);
    [os_lat(:,1), os_lat(:,2)] = ...
        algebra.rotateY(os_lat(:,1), os_lat(:,2), pitchAngle);
end