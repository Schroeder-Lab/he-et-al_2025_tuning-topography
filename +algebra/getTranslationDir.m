function [vecDir, dir_rad, dir_vis] = getTranslationDir(T, r)

% INPUTS
% T     [azimuth elevation], translation vector pointing into direction the
%       surrounding is moving during translation (diametrically to 
%       translation of subject/eye)
% r     [azimuth elevation], position of receptive field or point on
%       unitsphere at which translational flow is measured

% transform translation vector to cartesian coordinates
% (dimensions: AP, ML, DV)
T = deg2rad(T);
T = [cos(T(2)) * cos(T(1)); ...
    cos(T(2)) * sin(T(1)); ...
    sin(T(2))];

% transform vector pointing to RF to cartesian coordinates
% (dimensions: AP, ML, DV)
r = deg2rad(r);
r_cart = [cos(r(2)) * cos(r(1)); ...
    cos(r(2)) * sin(r(1)); ...
    sin(r(2))];

% instantaneous optic flow at point r given by component of T that is
% tangential to the sphere
v = T - (T' * r_cart) .* r_cart;

% project flow vector onto tangent plane of r with sphere, use basis
% vectors parallel to horizontal and coronal plane (increasing elevation
% given by derivative with respect to elevation)
e_phi = [-sin(r(1)); cos(r(1)); 0]; % "horizontal" direction
e_theta = [-sin(r(2))*cos(r(1)); -sin(r(2))*sin(r(1)); cos(r(2))]; % "vertical" direction
vecDir = [v' * e_phi; v' * e_theta]; % projected flow vector

% determine flow direction
dir_rad = atan(vecDir(2) / vecDir(1));
if vecDir(1) < 0
    dir_rad = dir_rad + pi;
end
% transform to our convention: 0 deg - towards right (as usual), 90 deg -
% downwards (usually upward)
dir_rad = mod(2*pi - dir_rad, 2*pi);
% transform to degrees
dir_vis = rad2deg(dir_rad);


% % make plots
% figure
% % sphere
% [X,Y,Z] = sphere;
% surf(X, Y, Z, 'FaceAlpha', 0.5)
% axis equal
% set(gca, 'YDir', 'reverse')
% xlabel('X (AP)')
% ylabel('Y (ML)')
% zlabel('Z (DV)')
% 
% hold on
% % translation vector from origin
% l = 1.5;
% plot3(l*[0 T(1)], l*[0 T(2)], l*[0 T(3)], 'b', 'LineWidth', 2)
% 
% % vector to RF
% plot3([0 r_cart(1)], [0 r_cart(2)], [0 r_cart(3)], ...
%     'k.-', 'LineWidth', 2, 'MarkerSize', 30)
% 
% % 3D flow vector
% plot3(r_cart(1) + [0 v(1)], r_cart(2) + [0 v(2)], r_cart(3) + [0 v(3)], ...
%     'r.-', 'LineWidth', 2, 'MarkerSize', 20)
% 
% % vectors spanning tangent plane
% plot3(r_cart(1) + [0 e_phi(1)], r_cart(2) + [0 e_phi(2)], r_cart(3) + [0 0], ...
%     'r', 'LineWidth', 2)
% plot3(r_cart(1) + [0 e_theta(1)], r_cart(2) + [0 e_theta(2)], ...
%     r_cart(3) + [0 e_theta(3)], ...
%     'r', 'LineWidth', 2)
% 
% % plot flow vector projected into tengant plane
% figure
% hold on
% plot([-1 1], [0 0], 'k')
% plot([0 0], [-1 1], 'k')
% plot([0 vecDir(1)], [0 vecDir(2)], 'r', 'LineWidth', 2)
% axis equal tight