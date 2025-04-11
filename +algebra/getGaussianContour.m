function [x, y] = getGaussianContour(var1, var2)
% From: https://waterprogramming.wordpress.com/2016/11/07/plotting-probability-ellipses-for-bivariate-normal-distributions/

valid = ~any(isnan([var1 var2]), 2);
var1 = var1(valid);
var2 = var2(valid);

% means
m1 = mean(var1);
m2 = mean(var2);
% covariance
CV = cov(var1, var2);
% eigen vectors and values
[eiVec, eiVal] = eig(CV);

% angles
theta = 0 : 0.01 : 2*pi;

% rotation matrix
% x_vec = [1; 0];
% cosRot = dot(x_vec, eiVec(:,1)) / (norm(x_vec) * norm(eiVec(:,1)));
rotation = pi/2 - acos(eiVec(1,1));
R = [sin(rotation) cos(rotation); -cos(rotation) sin(rotation)];

% ellipse
xRadius = 2 * sqrt(eiVal(1,1));
yRadius = 2 * sqrt(eiVal(2,2));
lineX = xRadius .* cos(theta);
lineY = yRadius .* sin(theta);
rotated = R * [lineX; lineY];
x = rotated(1,:)' + m1;
y = rotated(2,:)' + m2;

% figure
% scatter(var1, var2, 'k', "filled")
% hold on
% for k = 1:2
%     plot(m1 + [0 eiVec(k,1).*eiVal(k,k)], m2 + [0 eiVec(k,2).*eiVal(k,k)], ...
%         'LineWidth', 2)
% end
% axis equal
% plot(x, y, 'r', 'LineWidth', 2)