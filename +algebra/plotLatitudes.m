function plotLatitudes(axisPole, numLat, colour)

% axisPole  [azimuth, elevation], in degrees
% numLat
% colour

nPoints = 200;

% Convert to radians
phi_A = deg2rad(axisPole(1));
theta_A = deg2rad(axisPole(2));

% Compute axis unit vector A
A = [cos(theta_A) * cos(phi_A);
     cos(theta_A) * sin(phi_A);
     sin(theta_A)];

% Find two orthonormal vectors perpendicular to A (for defining plane 
% parallel to latitudes)
% One arbitrary orthogonal vector (not parallel to A)
if abs(dot(A, [1;0;0])) < 0.9
    temp = [1; 0; 0];  % use x-axis unless A is close to AP direction
else
    temp = [0; 0; 1];
end

v1 = cross(A, temp);
v1 = v1 / norm(v1);    % unit vector orthogonal to A
v2 = cross(A, v1);     % second vector orthogonal to both A and v1

% Generate points along multiple longitudes
theta = linspace(0, 2 * pi, nPoints);  % parameter around circle

axis equal;
xlabel('X'); ylabel('Y'); zlabel('Z');
view(3);
grid on;
set(gca, 'YDir', 'reverse')

% Plot the axis
quiver3(0,0,0, A(1), A(2), A(3), 1.2, 'LineWidth', 2, 'Color', colour);
quiver3(0,0,0, -A(1), -A(2), -A(3), 1.2, 'LineWidth', 2, 'Color', colour);

% Loop through and plot longitudes
for k = 1:numLat
    % angle from pole (0 = pole, pi/2 = equator relative to A)
    angle = pi * k / (numLat + 1);

    % radius of latitude circle
    r = sin(angle);

    % centre of circle
    centre = cos(angle) * A;

    % Generate circle points in the plane perpendicular to A
    points = zeros(3, nPoints);
    for i = 1:nPoints
        points(:,i) = centre + r * ...
            (cos(theta(i)) * v1 + sin(theta(i)) * v2);
    end

    % Plot latitude
    plot3(points(1,:), points(2,:), points(3,:), 'Color', colour, 'LineWidth', 1);
end
