function [vecOri, ori_rad, ori_vis] = getLatitudeOrientation(A, r)

% INPUTS
% A     [azimuth elevation], vector which is perpendicular to latitudes
% r     [azimuth elevation], position of receptive field or point on
%       unitsphere at which translational flow is measured

% transform vector spanning latitudes to cartesian coordinates
% (dimensions: AP, ML, DV)
A = deg2rad(A);
A = [cos(A(2)) * cos(A(1)); ...
    cos(A(2)) * sin(A(1)); ...
    sin(A(2))];

% transform vector pointing to RF to cartesian coordinates
% (dimensions: AP, ML, DV)
r = deg2rad(r);
r_cart = [cos(r(2)) * cos(r(1)); ...
    cos(r(2)) * sin(r(1)); ...
    sin(r(2))];

% orientation of latitude circle
v = cross(r_cart, A);

% project orientation vector onto tangent plane of r with sphere, use basis
% vectors parallel to horizontal and coronal plane (increasing elevation
% given by derivative with respect to elevation)
e_phi = [-sin(r(1)); cos(r(1)); 0]; % "horizontal" direction
e_theta = [-sin(r(2))*cos(r(1)); -sin(r(2))*sin(r(1)); cos(r(2))]; % "vertical" direction
vecOri = [v' * e_phi; v' * e_theta]; % projected orientation vector

% determine orientation
ori_rad = atan(vecOri(2) / vecOri(1));
% transform to our convention
ori_rad = mod(2*pi - ori_rad, pi);
% transform to degrees
ori_vis = rad2deg(ori_rad);


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
% plot3(l*[0 A(1)], l*[0 A(2)], l*[0 A(3)], 'b', 'LineWidth', 2)
% 
% % vector to RF
% plot3([0 r_cart(1)], [0 r_cart(2)], [0 r_cart(3)], ...
%     'k.-', 'LineWidth', 2, 'MarkerSize', 30)
% 
% % 3D orientation vector
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
% plot([0 vecOri(1)], [0 vecOri(2)], 'r', 'LineWidth', 2)
% axis equal tight