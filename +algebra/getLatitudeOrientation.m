function [vecOri, ori_rad, ori_vis] = getLatitudeOrientation(A, r)
%GETLATITUDEORIENTATION   Determine preferred orientation predicted from
%latitudes at RF location.

% INPUTS
% A         [azimuth elevation], vector which is perpendicular to latitudes
% r         [azimuth elevation], position of receptive field or point on
%           unitsphere at which translational flow is measured

% OUTPUTS
% vecOri    [x y], orientation vector projected onto tangent plane of RF 
%           with sphere
% ori_rad   double, preferred orientation in radian
% ori_vis   double, preferred orientation in degrees

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