function [azimuth2, elevation2] = rotateY(azimuth, elevation, alpha)
%ROTATEY   Determine new vector position on sphere after changing pitch of
%sphere (rotation around y-axis).

% INPUTS
% azimuth       double, azimuth of vector
% elevation     double, elevation of vector
% alpha         double, rotation angle in degrees

% OUTPUTS
% azimuth2      double, azimuth of vector after rotation
% elevation2    double, elevation of vector after rotation

% Convert to radians
azimuth = deg2rad(azimuth);
elevation = deg2rad(elevation);
alpha = deg2rad(alpha);

% Cartesian coordinates on unit sphere
x = cos(elevation) .* cos(azimuth);
y = cos(elevation) .* sin(azimuth);
z = sin(elevation);

% Rotation around +y axis (right-hand rule)
x2 =  x .* cos(alpha) + z .* sin(alpha);
y2 =  y;
z2 = -x .* sin(alpha) + z .* cos(alpha);

% Back to azimuth/elevation
azimuth2 = atan2(y2, x2);
elevation2 = asin(z2);

% Convert back to degrees
azimuth2 = rad2deg(azimuth2);
elevation2 = rad2deg(elevation2);