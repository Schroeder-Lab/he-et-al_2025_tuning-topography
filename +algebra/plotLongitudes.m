function plotLongitudes(axisPole, numLong, colour)
%PLOTLONGITUDES   Plot longitudinal fields on sphere given specific axis.

% axisPole  [azimuth, elevation], axis in degrees
% numLat    integer, number of longitudes to plot
% colour    colour of longitudes

nPoints = 200;

% Convert to radians
phi_A = deg2rad(axisPole(1));
theta_A = deg2rad(axisPole(2));

% Compute axis unit vector A
A = [cos(theta_A) * cos(phi_A);
     cos(theta_A) * sin(phi_A);
     sin(theta_A)];

% Find two orthonormal vectors perpendicular to A (to sweep out the longitudes)
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
theta = linspace(0, pi, nPoints);  % parameter around semi-circles

% figure;
% hold on;
axis equal;
xlabel('X'); ylabel('Y'); zlabel('Z');
view(3);
grid on;
set(gca, 'YDir', 'reverse')

% Plot the axis
quiver3(0,0,0, A(1), A(2), A(3), 1.2, 'LineWidth', 2, 'Color', colour);
quiver3(0,0,0, -A(1), -A(2), -A(3), 1.2, 'LineWidth', 2, 'Color', colour);

% Loop through and plot longitudes
for k = 0:(numLong-1)
    angle = 2*pi * k / numLong;

    % plane defined by A and rotated v1
    dir = cos(angle)*v1 + sin(angle)*v2;

    % generate circle points: rotate dir around A
    points = zeros(3, nPoints);
    for i = 1:nPoints
        t = theta(i);
        points(:,i) = cos(t)*A + sin(t)*dir;
    end

    % Plot longitude
    plot3(points(1,:), points(2,:), points(3,:), 'Color', colour, 'LineWidth', 1);
end
