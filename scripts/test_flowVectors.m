% receptive field position
phi_deg = 80; % azimuth
theta_deg = 20; % elevation
% translation vector
T = [-1; 0; 0]; % forward translation

% transform phi to angle from 90 deg azimuth (where 0 deg usually is)
% phi = 90 - phi_deg;

% convert degree to radian
phi = deg2rad(phi_deg);
theta = deg2rad(theta_deg);

% vector pointing to RF in Cartesian coordinate system 
% (dimensions: AP, ML, DV)
r = [cos(theta) * cos(phi); ...
    cos(theta) * sin(phi); ...
    sin(theta)];

% instantaneous optic flow at point r given by component of T that is
% tangential to the sphere
v = T - (T' * r) .* r;

% project flow vector onto tangent plane of r with sphere, use basis
% vectors parallel to horizontal and coronal plane (increasing elevation
% given by derivative with respect to elevation)
e_phi = [-sin(phi); cos(phi); 0]; % "horizontal" direction
e_theta = [-sin(theta)*cos(phi); -sin(theta)*sin(phi); cos(theta)]; % "vertical" direction
v_proj = [v' * e_phi; v' * e_theta]; % projected flow vector

% determine flow direction
direction = rad2deg(atan(v_proj(2) / v_proj(1)));
% transform to our convention
direction = 360 - direction;


% make plots
figure
% sphere
[X,Y,Z] = sphere;
surf(X, Y, Z, 'FaceAlpha', 0.5)
axis equal
set(gca, 'YDir', 'reverse')
xlabel('X (AP)')
ylabel('Y (ML)')
zlabel('Z (DV)')

hold on
% translation vector from origin
plot3(2*[0 T(1)], 2*[0 T(3)], 2*[0 T(3)], 'b', 'LineWidth', 2)

% vector to RF
plot3([0 r(1)], [0 r(2)], [0 r(3)], 'k.-', 'LineWidth', 2, 'MarkerSize', 30)

% 3D flow vector
plot3(r(1) + [0 v(1)], r(2) + [0 v(2)], r(3) + [0 v(3)], 'r.-', ...
    'LineWidth', 2, 'MarkerSize', 20)

% vectors spanning tangent plane
plot3(r(1) + [0 e_phi(1)], r(2) + [0 e_phi(2)], r(3) + [0 0], ...
    'r', 'LineWidth', 2)
plot3(r(1) + [0 e_theta(1)], r(2) + [0 e_theta(2)], r(3) + [0 e_theta(3)], ...
    'r', 'LineWidth', 2)

% plot flow vector projected into tengant plane
figure
hold on
plot([-1 1], [0 0], 'k')
plot([0 0], [-1 1], 'k')
plot([0 v_proj(1)], [0 v_proj(2)], 'r', 'LineWidth', 2)
axis equal tight
