function [azimuth2, elevation2] = rotateY(azimuth, elevation, alpha)

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